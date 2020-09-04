#include "AddOns/Analysis/Analyses/Analysis_Base.H"
#include "Analysis/Observable_Base.H"
#include "Analysis/ChannelAlgorithms/ChannelAlgorithm_Base.H"
#include "ATOOLS/Org/Gzip_Stream.H"
#include "Cumulant.H"

namespace RESUM {

  class NLO_Analysis_2: public ANALYSIS::Analysis_Base {
  public:

    NLO_Analysis_2(const ANALYSIS::Argument_Matrix &params);

    ~NLO_Analysis_2();

    double sumW = 0.;
    int n;
    void Evaluate(double weight, double ncount,int mode) override;
    void EndEvaluation(double scale) override;
    void EndEvaluation(std::vector<ATOOLS::Histogram*>& hists, double scale) override;
    void Output(const std::string& pname) override;
    void EvaluateNLOevt() override;

    void Fill(int i, double value, double weight, double ncount, int mode);
    std::vector<std::string> m_channels;
    Primitive_Observable_Base * Copy() const;
    std::string m_addition = "";
    std::string m_fname;
    int m_nborn = -1;
    int m_fills = 0;
    int m_fills_tmp = 0;
    bool m_recycle = true;
    bool m_oneFile = false;
    bool m_properZeroFill = false;
    bool m_filledOnce = false;

    void PrintHistos() {
      for(auto& h: m_histos) {
        msg_Out()<<h->Name()<<" has bins:"<<"\n";
        for(int i=0; i<h->Nbin(); i++)
          msg_Out()<<h->LowEdge(i)<<" "<<h->HighEdge(i)<<" "<<h->Bin(i)<<"\n";
      }
      
      msg_Out()<<"\n";
    }

    void AddZeroWeight(int ncount, int mode) {
      for(auto s: m_sigmas) {
        if(mode==0) s->AddZeroPoint(ncount);
        else if(mode==1) s->AddZeroPointMCB(ncount);
      }
      AddZeroPoint(ncount,mode);
    }

    
  private:

    ANALYSIS::Argument_Matrix m_params;

    std::vector<Observable_Base*> m_obss;

    std::vector<Cumulant::Ptr> m_sigmas;


    // std::map<std::string,std::vector<std::vector<double>>> m_sigma;
    
    double EtaBeam(const ATOOLS::Vec4D& p, size_t beamId);
    
    std::vector<ChannelAlgorithm_Base::Ptr> m_channelAlgs;

    void Output(My_Out_File& ofile);
    void Output(ATOOLS::Gzip_Stream& ogzip);
  };// end of class NLO_Analysis_2

}// end of namespace RESUM

#include "AddOns/Analysis/Main/Analysis_Handler.H"
#include "SHERPA/Single_Events/Event_Handler.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Vector.H"
#include "Tools/StringTools.H"
#include <algorithm>


using namespace RESUM;
using namespace ANALYSIS;
using namespace ATOOLS;

NLO_Analysis_2::NLO_Analysis_2(const Argument_Matrix &params):
  Analysis_Base(params[0][0]), m_params(params)
{
  DEBUG_FUNC(this);
  m_name+="_Resum";
  Data_Reader reader(",",";","!","=");
  Algebra_Interpreter *ip=reader.Interpreter();
  m_channelAlgs.clear();

  // msg_Out()<<"Start\n";
  // ATOOLS::Blob* sig = Analysis()->AnalysisHandler()->EventHandler()->GetBlobs()->FindFirst(btp::Signal_Process);
  // if(not sig) msg_Out()<<"No signal.\n";
  // else msg_Out()<<"found sig\n";
  // auto* vars = &(*sig)["Variation_Weights"]->Get<SHERPA::Variation_Weights>();
  // if(not vars) msg_Out()<<"No vars.\n";
  // const size_t nvars = vars->GetNumberOfVariations();
  // std::map<std::string,double> varweights;
  // if(nvars != 0) {
  //   for(size_t i=0; i<nvars; i++) {
  //     varweights[vars->GetVariationNameAt(i)] = vars->GetVariationWeightAt(i, SHERPA::Variations_Type::main);
  //   }
  //   SetVarWeights(&varweights);
  // }

  for (size_t i(1);i<params.size();++i) {
    std::vector<std::string> ps = {params[i].begin()+1,
                                   params[i].end()};
    if (params[i][0] == "Options") {
      Key_Base opts("Options",params[i]);
      m_recycle = RESUM::to_type<bool>(opts.KwArg("REUSE_ALGS","1"));
      m_oneFile = RESUM::to_type<bool>(opts.KwArg("ONEFILE","0"));
      if(m_oneFile) m_fname = opts.KwArg("FILENAME","ResumResults");
      continue;
    }
    if (params[i][0] == "ChAlg") {
      m_channelAlgs.emplace_back(ChAlg_Getter::GetObject
                                 (params[i][1], {params[i][1],ps}));
      for(const std::string& ch: m_channelAlgs.back()->ChannelNames(true))
        m_channels.push_back(ch);
      continue;
    }
  }
  for (size_t i(1);i<params.size();++i) {
    if (params[i][0] == "Options") continue;
    if (params[i][0] == "ChAlg") continue;
    if (params[i].size()<5) continue;    
    std::vector<std::string> ps = {params[i].begin()+1,
                                   params[i].end()};
    Observable_Base *obs(Observable_Getter::GetObject
			 (params[i][0],{params[i][0], ps}));
    if (obs == nullptr) {
      msg_Error()<<METHOD<<"(): Observable not found '"<<params[i][0]<<"'.\n";
      continue;
    }
    std::valarray<double> edges;
    size_t pos = params[i][1].find("EDGES:");
    if (pos != std::string::npos) {
      const std::vector<std::string> es = RESUM::split(params[i][1].substr(pos+6),"_");
      edges = std::valarray<double>(es.size());
      for(size_t j=0; j<es.size(); j++) edges[j] = ToType<double>(ip->Interprete(es[j]));
      
    }
    else {
      const double xmin = ToType<double>(ip->Interprete(params[i][1]));
      const double xmax = ToType<double>(ip->Interprete(params[i][2]));
      const size_t nbin = ToType<size_t>(ip->Interprete(params[i][3]));
      const int tp = HistogramType(params[i][4]);
      edges = Cumulant::Edges(tp,xmin,xmax,nbin);
      msg_Debugging()<<"Init '"<<params[i][0]<<"', type "<<tp<<" ("<<params[i][4]<<") "
                     <<" with "<<nbin<<" bins in ["<<xmin<<","<<xmax<<"]\n";
    }
    m_sigmas.push_back(std::make_shared<Cumulant>(edges, obs->Tag()+"_Sigma"));
    // m_sigmas.back()->InitVarWeights(p_varweights);
    // m_sigmas.back()->SetVariants(m_channels);
    m_obss.push_back(obs);
  }
}

NLO_Analysis_2::~NLO_Analysis_2()
{
  for (size_t i(0);i<m_obss.size();++i) delete m_obss[i];
}

double NLO_Analysis_2::EtaBeam(const Vec4D& p, size_t beamId) {
  double eta = 1./(2.0*rpa->gen.PBeam(beamId)[0]);
  if(beamId == 0) {
    eta *= p.PPlus();
  }
  else if(beamId == 1) {
    eta *= p.PMinus();
  }
  else THROW(fatal_error, "Something went wrong, beamId > 1.");
  return eta;
}

void NLO_Analysis_2::Evaluate(double weight,double ncount,int mode)
{
  DEBUG_FUNC("mode = "<<mode);
  // if(Analysis()->Sub() and Analysis()->Real() and
  //    Analysis()->Sub() != Analysis()->Real()) {
  //   msg_Debugging()<<"Subtraction.\n";
  // }
  // else {
  //   msg_Debugging()<<"Real.\n";
  // }
  // msg_Debugging()<<"Weight = "<<weight<<" "<<mode<<"\n";
  // if(p_varweights) {
  //   for(const auto& v: *p_varweights) msg_Debugging()<<v.first<<" "<<v.second<<"\n";
  // }
  // else {
  //   msg_Debugging()<<"No variation weights.\n";
  // }
  // msg_Debugging()<<"\n";
  // sumW += weight;
  // n += ncount;
  // msg_Out()<<"weight: "<<sumW/n<<"\n";
  // calculate observable for event and fill into histo
  if(!m_filledOnce) {
    for(auto& s: m_sigmas) {
      if(p_varweights) s->InitVarWeights(p_varweights);
      s->SetVariants(m_channels);
    }
    m_filledOnce=true;
  }
  Particle_List all(*p_ana->GetParticleList(m_listname));
  Vec4D_Vector mom(2+all.size());
  Flavour_Vector fl(2+all.size());
  for (size_t i=0; i<all.size(); i++) {
    mom[2+i]=all[i]->Momentum();
    fl[2+i]=all[i]->Flav();
  }
  if(Analysis()->Sub() and Analysis()->Real() and
     Analysis()->Sub() != Analysis()->Real()) {
    // for S-Event get mapped flavours
    // @TODO: could also do it this way for R, what is the proper way
    fl[0] = Analysis()->Sub()->p_fl[0];
    fl[1] = Analysis()->Sub()->p_fl[1];
  }
  else {
    // R, VI and B envents
    const ATOOLS::Blob* sig = Analysis()->AnalysisHandler()->EventHandler()->GetBlobs()->FindFirst(btp::Signal_Process);
    fl[0] = sig->ConstInParticle(0)->Flav();
    fl[1] = sig->ConstInParticle(1)->Flav();
  }

  std::vector<std::string> ch(m_channelAlgs.size(),"");
  for(size_t i=0; i<ch.size(); i++) {
    ch[i] = m_channelAlgs[i]->Channel(mom,fl,2,true);
    if(std::find(m_channels.begin(), m_channels.end(), ch[i]) == m_channels.end()) {
      msg_Error()<<"Channel not found: "<<ch[i]<<"\n";
      msg_Error()<<"Available:\n";
      for(auto& c: m_channels) msg_Error()<<c<<"\n";
    } 
  }


  std::map<std::string, typename Algorithm<double>::Ptr> algorithms;
  for (size_t i=0; i<m_obss.size(); i++) {
    const double value = m_recycle ?  m_obss[i]->Value(mom, fl, algorithms) : m_obss[i]->Value(mom, fl);
    msg_Debugging()<<"value["<<i<<"] = "<<value
                   <<" ( w = "<<weight<<", n = "<<ncount<<" )\n";
    if(mode==0) {
      m_sigmas[i]->Fill(value,weight,p_varweights,ncount);
      m_sigmas[i]->Fill(ch);
    }
    else if(mode==1) {
      m_sigmas[i]->FillMCB(value,weight,p_varweights,ncount);
      m_sigmas[i]->FillMCB(ch);
    }
    else THROW(fatal_error,"Unknown mode.");
  }
  return;
}

void NLO_Analysis_2::Fill(int i, double value, double weight, double ncount, int mode) {
  if(mode==0) m_sigmas[i]->Fill(value,weight,p_varweights,ncount);
  else if(mode==1) m_sigmas[i]->FillMCB(value,weight,p_varweights,ncount);
  else THROW(fatal_error,"Unknown mode.");
}

void NLO_Analysis_2::EvaluateNLOevt() {
  for(auto& s: m_sigmas) s->FinishMCB();
  Analysis_Base::EvaluateNLOevt();
}



void NLO_Analysis_2::EndEvaluation(std::vector<ATOOLS::Histogram*>& hists, double scale)
{
  // for(auto& h: hists) {
  //   h->MPISync();
  // }
  // auto histos = hists;
  // hists.clear();
  // for (size_t i=0; i<histos.size(); i+=2+histos[i]->Nbin()) {
  //   hists.push_back(histos[i]);
  //   m_sigma[histos[i]->Name()] = std::vector<std::vector<double>>(histos[i]->Nbin()+1);


  //   auto& h = histos[i+1];
  //   if(h->Nbin()!=2) THROW(fatal_error,"All histograms here should have two bins, this has "+std::to_string(h->Nbin()));

  //   double v = histos[i]->LowEdge(0);
  //   double N = m_properZeroFill ? h->Fills() : m_fills;
  //   double sigW = h->Bin(1)/N;
  //   double sigW2 = h->Bin2(1)/N;
  //   double sigErr = sqrt((sigW2-sqr(sigW))/(N-1));
  //   // msg_Out()<<v<<" "<<sigW<<" "<<sigW2<<" "<<sigErr<<" ";
  //   double BsigW = h->Bin(2)/N;
  //   double BsigW2 = h->Bin2(2)/N;
  //   double BsigErr = sqrt((BsigW2-sqr(BsigW))/(N-1));
  //   // msg_Out()<<BsigW<<" "<<BsigW2<<" "<<BsigErr<<"\n";
  //   m_sigma[histos[i]->Name()][0] = {v, sigW, sigW2, sigErr, BsigW, BsigW2, BsigErr, N};
  //   for(int j=0; j<histos[i]->Nbin(); j++) {
  //     auto& h = histos[i+j+2];
  //     if(h->Nbin()!=2) THROW(fatal_error,"All histograms here should have two bins, this has "+std::to_string(h->Nbin()));
  //     v = histos[i]->HighEdge(j);
  //     N = m_properZeroFill ? h->Fills() : m_fills;
  //     sigW = h->Bin(1)/N;
  //     sigW2 = h->Bin2(1)/N;
  //     sigErr = sqrt((sigW2-sqr(sigW))/(N-1));
  //     // msg_Out()<<v<<" "<<sigW<<" "<<sigW2<<" "<<sigErr<<" ";
  //     BsigW = h->Bin(2)/N;
  //     BsigW2 = h->Bin2(2)/N;
  //     BsigErr = sqrt((BsigW2-sqr(BsigW))/(N-1));
  //     // msg_Out()<<BsigW<<" "<<BsigW2<<" "<<BsigErr<<"\n";
  //     m_sigma[histos[i]->Name()][j+1] = {v, sigW, sigW2, sigErr, BsigW, BsigW2, BsigErr, N};
  //     delete histos[i+j+1];
  //   }
  // }
  Analysis_Base::EndEvaluation(hists, scale);
}

void NLO_Analysis_2::EndEvaluation(double scale) {
  for(auto& s: m_sigmas) s->Finalize();
  Analysis_Base::EndEvaluation(scale);
}

void NLO_Analysis_2::Output(My_Out_File& ofile) {
  // for(auto& sig: m_sigma) {
  //   //    msg_Out()<<sig.first<<"\n";
  //   *ofile<<"# "<<m_addition<<sig.first<<"_Sigma.dat\n";
  //   *ofile<<"# v Sigma{sumW sumW2 err} barSigma{sumW sumW2 err} NumEntries\n";
  //   for(auto& vals: sig.second) {
  //     for(double val: vals) *ofile<<val<<" ";
  //     *ofile<<"\n";
  //   }
  //   *ofile<<"\n\n";
  // }
}

void NLO_Analysis_2::Output(Gzip_Stream& ofile) {
  for(auto& s: m_sigmas) {
    //    msg_Out()<<sig.first<<"\n";
    // this becomes title
    s->Output(ofile, "", "Channel");
    *ofile.stream()<<"\n\n";
  }
}

void NLO_Analysis_2::Output(const std::string& pname) {
// #ifdef USING__MPI
//   if (MPI::COMM_WORLD.Get_rank()) return;
// #endif
  if(m_oneFile) {
    std::string fname = pname+"/"+m_fname;
    //while(fname.back()=='/') fname.pop_back();
    //std::replace( fname.begin(), fname.end(), '/', '_');
    fname += ".dat";
    msg_Out()<<"Write to file "<<fname<<"\n";
    // My_Out_File ofile(fname);
    // ofile.Open();
    // ofile->precision(ToType<int>(rpa->gen.Variable("HISTOGRAM_OUTPUT_PRECISION")));
    // Output(ofile);
    // ofile.Close();
    fname += ".gz";
    Gzip_Stream ogzip;
    ogzip.open(fname);
    Output(ogzip);
    ogzip.close();
  }
  // else {
  //   for(auto& ch: m_channels) {
  //     ch.second->m_addition = "Channel_"+ch.first+"_";
  //     ch.second->Output(pname);
  //   }
  //   for(auto& sig: m_sigma) {
  //     My_Out_File ofile(pname+"/"+m_addition+sig.first+"_Sigma.dat");
  //     ofile.Open();
  //     ofile->precision(ToType<int>(rpa->gen.Variable("HISTOGRAM_OUTPUT_PRECISION")));
  //     *ofile<<"# v Sigma{sumW sumW2 err} barSigma{sumW sumW2 err} NumEntries\n";
  //     for(auto& vals: sig.second) {
  //       for(double val: vals) *ofile<<val<<" ";
  //       *ofile<<"\n";
  //     }
  //     ofile.Close();  
  //   }
  // }
  //Analysis_Base::Output(pname);
}

Primitive_Observable_Base *NLO_Analysis_2::Copy() const 
{
  return new NLO_Analysis_2(m_params);
}

DECLARE_GETTER(NLO_Analysis_2,"NLOResum2",
	       Primitive_Observable_Base,Argument_Matrix);

Primitive_Observable_Base *ATOOLS::Getter
<Primitive_Observable_Base,Argument_Matrix,NLO_Analysis_2>::
operator()(const Argument_Matrix &parameters) const
{
  if (parameters.size()==0 || parameters[0].size()<1) return nullptr;
  return new NLO_Analysis_2(parameters);
}

void ATOOLS::Getter
<Primitive_Observable_Base,Argument_Matrix,NLO_Analysis_2>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"list"; 
}
