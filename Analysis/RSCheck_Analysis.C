#include "AddOns/Analysis/Analyses/Analysis_Base.H"
#include "Analysis/Observable_Base.H"
#include "Analysis/ChannelAlgorithms/ChannelAlgorithm_Base.H"

namespace RESUM {

  class RSCheck_Analysis: public ANALYSIS::Analysis_Base {
  public:

    RSCheck_Analysis(const ANALYSIS::Argument_Matrix &params);

    ~RSCheck_Analysis();

    double sumW = 0.;
    int n;
    void Evaluate(double weight, double ncount,int mode) override;
    void EndEvaluation(double scale) override;
    void Output(const std::string& pname) override;
    void EvaluateNLOevt() override;

    void Fill(int i, double value, double weight, double ncount, int mode);
    std::map<std::string, RSCheck_Analysis*> m_channels;
    Primitive_Observable_Base * Copy() const;
    std::string m_addition = "";
    int m_nborn = -1;
    int m_fills = 0;
    int m_fills_tmp = 0;


    void PrintHistos() {
      for(auto& h: m_histos) {
        msg_Out()<<h->Name()<<" has bins:"<<"\n";
        for(int i=0; i<h->Nbin(); i++)
          msg_Out()<<h->LowEdge(i)<<" "<<h->HighEdge(i)<<" "<<h->Bin(i)<<"\n";
      }
      
      msg_Out()<<"\n";
    }
    
  private:

    ANALYSIS::Argument_Matrix m_params;

    std::vector<Observable_Base*> m_obss;

    std::map<std::string,std::vector<std::vector<double>>> m_sigma;
    
    double EtaBeam(const ATOOLS::Vec4D& p, size_t beamId);
    
    std::vector<ChannelAlgorithm_Base::Ptr> m_channelAlgs;
  };// end of class RSCheck_Analysis

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

RSCheck_Analysis::RSCheck_Analysis(const Argument_Matrix &params):
  Analysis_Base(params[0][0]), m_params(params)
{
  DEBUG_FUNC(this);
  m_name+="_Resum";
  Data_Reader reader(",",";","!","=");
  Algebra_Interpreter *ip=reader.Interpreter();
  m_channelAlgs.clear();
  for (size_t i(1);i<params.size();++i) {
    std::vector<std::string> ps = {params[i].begin()+1,
                                   params[i].end()};
    if (params[i][0] == "ChAlg") {
      m_channelAlgs.emplace_back(ChAlg_Getter::GetObject
                                 (params[i][1], {params[i][1],ps}));
      continue;
    }
    if (params[i].size()<5) continue;    
    Observable_Base *obs(Observable_Getter::GetObject
			 (params[i][0],{params[i][0], ps}));
    if (obs == nullptr) {
      msg_Error()<<METHOD<<"(): Observable not found '"<<params[i][0]<<"'.\n";
      continue;
    }
    double xmin = ToType<double>(ip->Interprete(params[i][1]));
    double xmax = ToType<double>(ip->Interprete(params[i][2]));
    size_t nbin = ToType<size_t>(ip->Interprete(params[i][3]));
    int tp = HistogramType(params[i][4]);
    msg_Debugging()<<"Init '"<<params[i][0]<<"', type "<<tp<<" ("<<params[i][4]<<") "
		   <<" with "<<nbin<<" bins in ["<<xmin<<","<<xmax<<"]\n";
    m_histos.push_back(new Histogram(tp,xmin,xmax,nbin,obs->Tag()));
    m_histos.push_back(new Histogram(HistogramType("LinErr"),0,1,2,obs->Tag()+"_Sigma"));
    for(int j=0; j<nbin; j++)
      m_histos.push_back(new Histogram(HistogramType("LinErr"),0,1,2,obs->Tag()+"_Sigma"));
    m_obss.push_back(obs);
  }
  // to stay consistent with previous defaults
  if(m_channelAlgs.size()==0 &&
     reader.GetValue<std::string>("RESUM::SORT_CHANNELS","YES")!="NO") {
    int nborn = reader.GetValue<int>("RESUM::NBORN",-1);
    if(nborn>0) {
      m_channelAlgs.emplace_back(ChAlg_Getter::GetObject("KT2_ee",
                                                         {"KT2_ee",
                                                          {std::to_string(nborn),
                                                           ""}}));
      m_channelAlgs.emplace_back(ChAlg_Getter::GetObject("KT2_ee",
                                                         {"KT2_ee",
                                                          {std::to_string(nborn),
                                                           "BLAND"}}));
      m_channelAlgs.emplace_back(ChAlg_Getter::GetObject("KT2_ee",
                                                         {"KT2_ee",
                                                          {std::to_string(nborn),
                                                           "BLAND_Z"}}));
    }
  }
  // get channels
  if(rpa->gen.Variable("RESUM::ANALYSIS_DO_CHANNELS")!="NO") {
    rpa->gen.SetVariable("RESUM::ANALYSIS_DO_CHANNELS","NO");
    for(ChannelAlgorithm_Base::Ptr alg: m_channelAlgs)
      for(const std::string& ch: alg->ChannelNames(true))
        m_channels.emplace(ch,dynamic_cast<RSCheck_Analysis*>(Copy()));
    rpa->gen.SetVariable("RESUM::ANALYSIS_DO_CHANNELS","YES"); 
  }
}

RSCheck_Analysis::~RSCheck_Analysis()
{
  for (size_t i(0);i<m_obss.size();++i) delete m_obss[i];
  for (auto& ch: m_channels) if(ch.second) delete ch.second;
}

double RSCheck_Analysis::EtaBeam(const Vec4D& p, size_t beamId) {
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

void RSCheck_Analysis::Evaluate(double weight,double ncount,int mode)
{
  DEBUG_FUNC("mode = "<<mode);
  // sumW += weight;
  // n += ncount;
  // msg_Out()<<"weight: "<<sumW/n<<"\n";
  // calculate observable for event and fill into histo
  ATOOLS::Blob* sig = Analysis()->AnalysisHandler()->EventHandler()->GetBlobs()->FindFirst(btp::Signal_Process);
  const size_t nout = sig->GetOutParticles().size();
  Vec4D_Vector mom(2+nout);
  Flavour_Vector fl(2+nout);
  for (size_t i=0; i<nout; i++) {
    mom[2+i]=sig->ConstOutParticle(i)->Momentum();
    fl[2+i]=sig->ConstOutParticle(i)->Flav();
  }
  fl[0] = sig->ConstInParticle(0)->Flav();
  fl[1] = sig->ConstInParticle(1)->Flav();

  // msg_Out()<<"\n\n";
  // for(auto m: mom) msg_Out()<<m<<"\n";
  // msg_Out()<<"\n\n";

  std::vector<std::string> ch(m_channelAlgs.size(),"");
  for(size_t i=0; i<ch.size(); i++) {
    ch[i] = m_channelAlgs[i]->Channel(mom,fl,2,true);
    if(m_channels.find(ch[i]) == m_channels.end()) {
      msg_Error()<<"Channel not found: "<<ch[i]<<"\n";
      msg_Error()<<"Available:\n";
      for(auto& c: m_channels) msg_Error()<<c.first<<"\n";
    } 
  }
  for(auto& c: m_channels) {
    bool found = false;
    for(const std::string& cname: ch)
      if(c.first == cname) {
        found = true;
        break;
      }
    if(!found) c.second->AddZeroPoint(ncount,mode);                          
  }

  int obsId = 0;
  for (size_t i=0; i<m_histos.size(); i+=2+m_histos[i]->Nbin()) {
    if(mode == 0) m_fills += ncount;
    if(mode == 1) m_fills_tmp = ncount;
    const double value = m_obss[obsId]->Value(mom, fl);
    msg_Debugging()<<"value["<<i<<"] = "<<value
                   <<" ( w = "<<weight<<", n = "<<ncount<<" )\n";
    // if(value < 1e-10) {
    //   msg_Out()<<value<<"\n";
    //   msg_Out()<<"\n\n";
    //   for(auto m: mom) msg_Out()<<m<<"\n";
    //   msg_Out()<<"\n\n";
    //   exit(1);
    // }
    Fill(i,value,weight,ncount,mode);
    for(const std::string& c: ch) {
      msg_Debugging()<<"Filling "<<c<<"\n";
      m_channels.at(c)->Fill(i,value,weight,ncount,mode);
    }
    obsId++;
  }
  return;
}

void RSCheck_Analysis::Fill(int i, double value, double weight, double ncount, int mode) {
  FillHisto(i,value,weight,ncount,mode);
  double fill = 0.75;
  if(value < m_histos[i]->LowEdge(0))
    fill = 0.25;
  FillHisto(i+1,fill,weight,ncount,mode);
  for(int j=0; j<m_histos[i]->Nbin(); j++) {
    fill = 0.75;
    if(value < m_histos[i]->HighEdge(j))
      fill = 0.25;
    msg_Debugging()<<m_histos[i]->HighEdge(j)<<" "<<value<<" "<<fill<<"\n";
    FillHisto(i+j+2,fill,weight,ncount,mode);
  }
}

void RSCheck_Analysis::EvaluateNLOevt() {
  m_fills += m_fills_tmp;
  m_fills_tmp = 0;
  for(auto& ch: m_channels) {
    ch.second->EvaluateNLOevt();
  }
  Analysis_Base::EvaluateNLOevt();
}

void RSCheck_Analysis::EndEvaluation(double scale)
{
  for(auto& ch: m_channels) {
    //msg_Out()<<ch.first<<"\n";
    ch.second->EndEvaluation(scale);
  }
  for(auto& h: m_histos) {
    h->MPISync();
  }
  auto histos = m_histos;
  m_histos.clear();
  for (size_t i=0; i<histos.size(); i+=2+histos[i]->Nbin()) {
    m_histos.push_back(histos[i]);
    m_sigma[histos[i]->Name()] = std::vector<std::vector<double>>(histos[i]->Nbin()+1);


    auto& h = histos[i+1];
    if(h->Nbin()!=2) THROW(fatal_error,"All histograms here should have two bins, this has "+std::to_string(h->Nbin()));

    double v = histos[i]->LowEdge(0);
    double N = h->Fills();
    double sigW = h->Bin(1)/N;
    double sigW2 = h->Bin2(1)/N;
    double sigErr = sqrt((sigW2-sqr(sigW))/(N-1));
    // msg_Out()<<v<<" "<<sigW<<" "<<sigW2<<" "<<sigErr<<" ";
    double BsigW = h->Bin(2)/N;
    double BsigW2 = h->Bin2(2)/N;
    double BsigErr = sqrt((BsigW2-sqr(BsigW))/(N-1));
    // msg_Out()<<BsigW<<" "<<BsigW2<<" "<<BsigErr<<"\n";
    m_sigma[histos[i]->Name()][0] = {v, sigW, sigW2, sigErr, BsigW, BsigW2, BsigErr, N};

    for(int j=0; j<histos[i]->Nbin(); j++) {
      auto& h = histos[i+j+2];
      if(h->Nbin()!=2) THROW(fatal_error,"All histograms here should have two bins, this has "+std::to_string(h->Nbin()));
      double v = histos[i]->HighEdge(j);
      double N = h->Fills();
      double sigW = h->Bin(1)/N;
      double sigW2 = h->Bin2(1)/N;
      double sigErr = sqrt((sigW2-sqr(sigW))/(N-1));
      // msg_Out()<<v<<" "<<sigW<<" "<<sigW2<<" "<<sigErr<<" ";
      double BsigW = h->Bin(2)/N;
      double BsigW2 = h->Bin2(2)/N;
      double BsigErr = sqrt((BsigW2-sqr(BsigW))/(N-1));
      // msg_Out()<<BsigW<<" "<<BsigW2<<" "<<BsigErr<<"\n";
      m_sigma[histos[i]->Name()][j+1] = {v, sigW, sigW2, sigErr, BsigW, BsigW2, BsigErr, N};
      delete histos[i+j+1];
    }
  }
  Analysis_Base::EndEvaluation(scale);
}

void RSCheck_Analysis::Output(const std::string& pname) {
// #ifdef USING__MPI
//   if (MPI::COMM_WORLD.Get_rank()) return;
// #endif
  for(auto& ch: m_channels) {
    ch.second->m_addition = "Channel_"+ch.first+"_";
    ch.second->Output(pname);
  }
  for(auto& sig: m_sigma) {
    My_Out_File ofile(pname+"/"+m_addition+sig.first+"_Sigma.dat");
    ofile.Open();
    ofile->precision(ToType<int>(rpa->gen.Variable("HISTOGRAM_OUTPUT_PRECISION")));
    *ofile<<"# v Sigma{sumW sumW2 err} barSigma{sumW sumW2 err} NumEntries\n";
    for(auto& vals: sig.second) {
      for(double val: vals) *ofile<<val<<" ";
      *ofile<<"\n";
    }                         
  }
  Analysis_Base::Output(pname);
}

Primitive_Observable_Base *RSCheck_Analysis::Copy() const 
{
  return new RSCheck_Analysis(m_params);
}

DECLARE_GETTER(RSCheck_Analysis,"RSCheck",
	       Primitive_Observable_Base,Argument_Matrix);

Primitive_Observable_Base *ATOOLS::Getter
<Primitive_Observable_Base,Argument_Matrix,RSCheck_Analysis>::
operator()(const Argument_Matrix &parameters) const
{
  if (parameters.size()==0 || parameters[0].size()<1) return nullptr;
  return new RSCheck_Analysis(parameters);
}

void ATOOLS::Getter
<Primitive_Observable_Base,Argument_Matrix,RSCheck_Analysis>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"list"; 
}
