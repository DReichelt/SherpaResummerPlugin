#include "SHERPA/Tools/Userhook_Base.H"
#include "SHERPA/Main/Sherpa.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "SHERPA/Initialization/Initialization_Handler.H" 
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Process/Process_Group.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "Main/Resum.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "Tools/Key_Base.H"
#include "Tools/StringTools.H"
#include "ATOOLS/Org/Gzip_Stream.H"

namespace RESUM {
class ResumHook : public SHERPA::Userhook_Base {
private:
  SHERPA::Sherpa* p_sherpa;
  RESUM::Resum* p_resum;
public:
  ResumHook(const SHERPA::Userhook_Arguments args) : Userhook_Base("Resum") {
    p_sherpa = args.p_sherpa;
    p_resum = new RESUM::Resum((*p_sherpa->GetInitHandler()->GetISRHandlers())[PDF::isr::id::hard_process], p_sherpa->GetInitHandler()->GetModel());
    p_resum->SetVariationParameters(p_sherpa->GetInitHandler()->GetVariations());
    
    for(PHASIC::Process_Base* p: p_sherpa->GetInitHandler()->GetMatrixElementHandler()->AllProcesses()) {
      for(size_t i=0; i<p->Size(); i++) {
        p_resum->AddOrder((*p)[i]->MaxOrder(0));
      }
    }

    ATOOLS::Algebra_Interpreter* ip = p_sherpa->GetInitHandler()->DataReader()->Interpreter();
    std::vector<Key_Base> variations;
    for(const std::vector<std::string>& params: GetParameters()) {
      std::vector<std::string> ps = {params.begin()+1, params.end()};
      if(params[0] == "ChAlg") p_resum->AddChAlg(params[1],ps);
      if(params[0] == "Var") variations.emplace_back(params[1],ps);
    }
    for(const std::vector<std::string>& params: GetParameters()) {
      if(params[0] == "ChAlg" or params[0] == "Var") continue;
      std::vector<std::string> ps = {params.begin()+1, params.end()};
      std::valarray<double> edges;
      size_t pos = params[1].find("EDGES:");
      if (pos != std::string::npos) {
        const std::vector<std::string> es = RESUM::split(params[1].substr(pos+6),"_");
        edges = std::valarray<double>(es.size());
        for(size_t j=0; j<es.size(); j++) edges[j] = to_type<double>(ip->Interprete(es[j])); 
      }
      else {
        if (params.size()<5) continue;
        const double xmin = to_type<double>(ip->Interprete(params[1]));
        const double xmax = to_type<double>(ip->Interprete(params[2]));
        const int nbin = to_type<int>(ip->Interprete(params[3]));
        const std::string& htype = params[4];
        msg_Debugging()<<"Going to add "<<params[0]<<" xmin = "<<xmin
                       <<" xmax = "<<xmax<<" Nbin = "<<nbin<<" Type = "<<htype<<"\n";
        edges = Cumulant::Edges(htype,xmin,xmax,nbin);
      }
      Observable_Key ObsKey = Observable_Key(params[0],ps);
      const std::string tag = ObsKey.KwArg("tag",ObsKey.Name());
      p_resum->AddObservable(ObsKey, edges, tag, true, true, true);
      for(const Key_Base& var: variations) {
        ObsKey = Observable_Key(params[0],ps);
        ObsKey.SetKwArg(var);
        if(var.Name() != "") {
          ObsKey.SetKwArg("tag",tag+"_"+var.Name());
        }
        const std::string& title = var.Name() == "" ?  tag : tag+" "+var.Name();
        p_resum->AddObservable(ObsKey, edges, title, true, false, false);
      }
    }
  }

  ~ResumHook() {
    delete p_resum;
  }

  ATOOLS::Return_Value::code Run(ATOOLS::Blob_List* blobs, double &weight) override {
    ATOOLS::Blob* sig = blobs->FindFirst(ATOOLS::btp::Signal_Process);
    if(not sig) THROW(fatal_error, "No signal process.");
    Treat(sig);
    return ATOOLS::Return_Value::Nothing;
  }

  void Treat(ATOOLS::Blob* sig) {
    const int ntrials = (*sig)["Trials"]->Get<double>();
    const double weight = (*sig)["Weight"]->Get<double>();
    const SHERPA::Matrix_Element_Handler* const me = p_sherpa->GetInitHandler()->GetMatrixElementHandler();
    if(not me) THROW(fatal_error, "No ME handler.");
    PHASIC::Process_Base* const proc = me->Process(); 
    if(not proc) THROW(fatal_error, "No process.");
    p_resum->SetVariationWeights(&(*sig)["Variation_Weights"]->Get<SHERPA::Variation_Weights>());
    if(proc->Info().Has(PHASIC::nlo_type::rsub)) {
      if(proc->GetSubevtList()) {
        for(ATOOLS::NLO_subevt* sevt: *proc->GetSubevtList()) {
          if(sevt->m_result == 0.) continue;
          ATOOLS::Cluster_Amplitude* ampl = ATOOLS::Cluster_Amplitude::New();
          ampl->SetNIn(sig->NInP());
          ampl->SetProc(proc);
          ampl->SetProcs(proc->AllProcs());
          for(size_t i=0; i<2; i++) {
            ampl->CreateLeg(-sevt->p_mom[i],sevt->p_fl[i].Bar());
          }
          for(size_t i=2; i<sevt->m_n; i++) {
            ampl->CreateLeg(sevt->p_mom[i],sevt->p_fl[i]);
          }
          ampl->SetMuF2(sevt->m_mu2[ATOOLS::stp::fac]);
          ampl->SetMuR2(sevt->m_mu2[ATOOLS::stp::ren]);
          ampl->SetMuQ2(sevt->m_mu2[ATOOLS::stp::res]);
          ampl->SetOrderQCD(proc->MaxOrder(0));
          ampl->SetOrderEW(proc->MaxOrder(1));
          p_resum->SetVariationIdx(sevt->m_idx);
          Treat(ampl,sevt->m_result);
          ampl->Delete();
        }
      }
    }
    else {
      ATOOLS::Cluster_Amplitude* ampl = ATOOLS::Cluster_Amplitude::New();
      for(const auto& p: sig->GetInParticles()) {
        ampl->CreateLeg(-p->Momentum(),p->Flav().Bar());
      }
      for(const auto& p: sig->GetOutParticles()) {
        ampl->CreateLeg(p->Momentum(),p->Flav());
      }
      ampl->SetNIn(sig->NInP());
      ampl->SetProc(proc);
      ampl->SetProcs(proc->AllProcs());
      ampl->SetMuF2(proc->ScaleSetter(true)->Scale(ATOOLS::stp::fac));
      ampl->SetMuR2(proc->ScaleSetter(true)->Scale(ATOOLS::stp::ren));
      ampl->SetMuQ2(proc->ScaleSetter(true)->Scale(ATOOLS::stp::res));
      ampl->SetOrderQCD(proc->MaxOrder(0));
      ampl->SetOrderEW(proc->MaxOrder(1));
      ampl->SetNLO(PHASIC::nlo_type::lo);
      p_resum->SetVariationIdx(-1);
      Treat(ampl,weight);
      ampl->Delete();
    }
    p_resum->FinishShowers(weight, ntrials);
  }

  void Treat(ATOOLS::Cluster_Amplitude* ampl, double weight) {
    p_resum->PrepareShower(ampl);
    p_resum->PerformShowers(weight);
    p_resum->CleanUp();
  }

  const std::vector<std::vector<std::string>>& GetParameters(bool strict=false) {
    if(strict or m_parameters.size() == 0) {
      std::string infile = p_sherpa->GetInitHandler()->File();
      if (infile.find('|')!=std::string::npos) infile=infile.substr(0,infile.find('|'));
      ATOOLS::Data_Reader localreader(" ",";","#");
      localreader.SetAddCommandLine(false);
      localreader.SetInputPath(p_sherpa->GetInitHandler()->Path()); 
      localreader.SetInputFile(infile+"|(resum){|}(resum)");
      localreader.MatrixFromFile(m_parameters,"");
    }
    return m_parameters;
  }


  void Finish() override {
    msg_Out()<<"Finish resummation.\n";
    std::string fname = "Results.dat.gz";
    m_ogzip.open(fname);
    for(auto c: p_resum->m_resultsNLL)     Finish(c.second, "NLL");
    for(auto c: p_resum->m_resultsExpLO)   Finish(c.second, "expLO");
    for(auto c: p_resum->m_resultsExpNLO)  Finish(c.second, "expNLO");
    for(auto c: p_resum->m_resultsLO)      Finish(c.second, "LO");
    for(auto c: p_resum->m_resultsNLO)     Finish(c.second, "NLO");
    m_ogzip.close();
  }


private:

  void Finish(Cumulant::Ptr cumulant, const std::string& title) {
    if(cumulant) {
      cumulant->Finalize();
      cumulant->Output(m_ogzip, title, title+" Channel");
      *m_ogzip.stream()<<"\n\n";
    }
  }


  std::vector<std::vector<std::string>> m_parameters = std::vector<std::vector<std::string>>(0);
  ATOOLS::Gzip_Stream m_ogzip;
};


}

using namespace SHERPA;
using namespace RESUM;

DECLARE_GETTER(ResumHook,"Resum",
               Userhook_Base,Userhook_Arguments);

Userhook_Base *ATOOLS::Getter<Userhook_Base,Userhook_Arguments,ResumHook>::
operator()(const Userhook_Arguments &args) const
{
  return new ResumHook(args);
}

void ATOOLS::Getter<Userhook_Base,Userhook_Arguments,ResumHook>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Resum userhook";
}

