#include "AddOns/Analysis/Analyses/Analysis_Base.H"

#include "Analysis/Observable_Base.H"
#include "Main/Resum.H"
#include "Analysis/ChannelAlgorithms/ChannelAlgorithm_Base.H"
#include "Analysis/Cumulant.H"
#include "ATOOLS/Org/Gzip_Stream.H"

namespace RESUM {

  class NLL_Analysis_2: public ANALYSIS::Analysis_Base {
  private:

    typedef std::vector<std::pair<double,double>> CumDist;
    typedef std::map<std::string,CumDist> ObsDist;
    typedef std::map<std::string,ObsDist> ChDist;
    
    ANALYSIS::Argument_Matrix m_params;

    std::vector<Observable_Base*> m_obss;

    Resum *p_resum;

    size_t m_nObs = 0;
    std::vector<std::string> m_ObsNames;
    std::map<std::string,Cumulant::Ptr> m_sigmas;
    double m_Ncount = 0;
    double m_sumW = 0;
    double m_sumW2 = 0;
    // ObsDist m_NLL;
    // ObsDist m_expLO;
    // ObsDist m_expNLO;
    // ChDist m_NLL_channels;
    // ChDist m_expLO_channels;
    // ChDist m_expNLO_channels;
    std::vector<std::string> m_channels;
    std::map<std::string,std::valarray<double>> m_expWeights;
    std::map<std::string,double> m_Weights;

    std::vector<ChannelAlgorithm_Base::Ptr> m_channelAlgs;    
  public:

    NLL_Analysis_2(const ANALYSIS::Argument_Matrix &params);

    void Evaluate(const ATOOLS::Blob_List & blobs,
		  double weight,double ncount) override;
    void Evaluate(double weight, double ncount,int mode) {}
    void Output(const std::string& pname) override;
    void EndEvaluation(double scale) override;
    Primitive_Observable_Base * Copy() const;

    void AddZeroWeight(int ncount, int mode) {
      for(auto s: m_sigmas) {
        if(mode==0) s.second->AddZeroPoint(ncount);
        else if(mode==1) s.second->AddZeroPointMCB(ncount);
      }
      AddZeroPoint(ncount,mode);
    }


  };// end of class NLL_Analysis_2

}// end of namespace RESUM

#include "AddOns/Analysis/Main/Analysis_Handler.H"
#include "SHERPA/Single_Events/Event_Handler.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"

#include <algorithm>

using namespace RESUM;
using namespace ANALYSIS;
using namespace ATOOLS;
using namespace std;

NLL_Analysis_2::NLL_Analysis_2(const Argument_Matrix& params):
  Analysis_Base(params[0][0]), m_params(params) {
  DEBUG_FUNC(this);
  m_name+="_Resum";
  ATOOLS::Data_Reader reader(",",";","!","=");
  ATOOLS::Algebra_Interpreter* ip = reader.Interpreter();
  p_resum=(RESUM::Resum*)ATOOLS::ToType<void*>(ATOOLS::rpa->gen.Variable("SHOWER_GENERATOR"));
  if (dynamic_cast<RESUM::Resum*>(p_resum) == nullptr) {
    THROW(fatal_error,"Resummer plugin not loaded");
  }
  p_resum->ResetObservables();
  m_channelAlgs.clear();
  m_Weights["LO"] = 0;
  m_Weights["NLO"] = 0;
  for(size_t i=1; i<params.size(); i++) {
    std::vector<std::string> ps = {params[i].begin()+1,
                                   params[i].end()};
    if (params[i][0] == "ChAlg") {
      m_channelAlgs.emplace_back(ChAlg_Getter::GetObject(params[i][1],
                                                         {params[i][1],
                                                          ps}));
      for(const std::string& ch: m_channelAlgs.back()->ChannelNames(true)) {
        m_channels.push_back(ch);
      }
      continue;
    }
  }
  for(size_t i=1; i<params.size(); i++) {
    if (params[i][0] == "ChAlg") continue;
    std::vector<std::string> ps = {params[i].begin()+1,
                                   params[i].end()};
    std::valarray<double> edges;
    size_t pos = params[i][1].find("EDGES:");
    if (pos != std::string::npos) {
      const std::vector<std::string> es = RESUM::split(params[i][1].substr(pos+6),"_");
      edges = std::valarray<double>(es.size());
      for(size_t j=0; j<es.size(); j++) edges[j] = ToType<double>(ip->Interprete(es[j])); 
    }
    else {
      if (params[i].size()<5) continue;
      const double xmin = ATOOLS::ToType<double>(ip->Interprete(params[i][1]));
      const double xmax = ATOOLS::ToType<double>(ip->Interprete(params[i][2]));
      const int nbin = ATOOLS::ToType<int>(ip->Interprete(params[i][3]));
      const string& htype = params[i][4];
      msg_Debugging()<<"Going to add "<<params[i][0]<<" xmin = "<<xmin
                     <<" xmax = "<<xmax<<" Nbin = "<<nbin<<" Type = "<<htype<<"\n";
      edges = Cumulant::Edges(htype,xmin,xmax,nbin);
    }
    // returns index of observable
    const string& obsName = p_resum->AddObservable({params[i][0],ps}, edges);
    m_nObs++;
    m_ObsNames.push_back(obsName);
    msg_Debugging()<<"Add "<<params[i][0]<<" as "<<obsName<<"\n";
    m_sigmas[obsName] = std::make_shared<Cumulant>(edges, obsName+"_Sigma", nullptr,1,"NLL");
    m_sigmas[obsName]->InitVarWeights(&m_Weights);
    m_sigmas[obsName]->SetVariants(m_channels);
  }
}

void NLL_Analysis_2::Evaluate(const ATOOLS::Blob_List& blobs,
			double weight, double ncount)
{
  DEBUG_FUNC("");
  if(!p_resum) THROW(fatal_error, "Plugin not loaded.");
  if(!p_resum->Initialized()) return;
  m_Ncount += ncount;
  m_sumW += weight;
  m_sumW2 += sqr(weight);
  Particle_List* allps = p_ana->GetParticleList(m_listname);
  if (allps == nullptr) AddZero(ncount,0);
  Particle_List all = *allps; 
  
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

  
  // find out channels
  vector<string> channels(m_channelAlgs.size());
  for(size_t i=0; i<channels.size(); i++) {
    channels[i] = m_channelAlgs[i]->Channel(mom,fl,2,true);
  }
  // get results from resum and fill
  m_Weights["NLL"] = weight;
  m_Weights["LO"] = weight;
  m_Weights["NLO"] = weight;
  for(size_t k=0; k<m_nObs; k++) {
    // msg_Out()<<k<<" "<<m_ObsNames.size()<<m_nObs<<"\n";
    const std::string& name = m_ObsNames[k];
    msg_Debugging()<<"Filling "<<name<<"\n";
    m_expWeights["NLL"] = p_resum->m_resNLL[k];
    m_expWeights["LO"] = p_resum->m_resExpLO[k];
    m_expWeights["NLO"] = p_resum->m_resExpNLO[k];
    m_sigmas[name]->Fill(&m_Weights,&m_expWeights,ncount);
    m_sigmas[name]->Fill(channels);
  }
}


void NLL_Analysis_2::EndEvaluation(double scale) {
  for(auto& s: m_sigmas) s.second->Finalize();
}

  
void NLL_Analysis_2::Output(const std::string& pname) {
  std::string fname = pname + "/" + "Results.dat";
  fname += ".gz";
  Gzip_Stream ogzip;
  ogzip.open(fname);
  for(auto& s: m_sigmas) {
    s.second->Output(ogzip, "", "Channel");
    *ogzip.stream()<<"\n\n";
  }
  ogzip.close();
}


Primitive_Observable_Base *NLL_Analysis_2::Copy() const 
{
  return new NLL_Analysis_2(m_params);
}

DECLARE_GETTER(NLL_Analysis_2,"NLLResum2",Primitive_Observable_Base,Argument_Matrix);

Primitive_Observable_Base *ATOOLS::Getter
<Primitive_Observable_Base,Argument_Matrix,NLL_Analysis_2>::
operator()(const Argument_Matrix &parameters) const
{
  if (parameters.size()==0 || parameters[0].size()<1) return nullptr;
  return new NLL_Analysis_2(parameters);
}

void ATOOLS::Getter
<Primitive_Observable_Base,Argument_Matrix,NLL_Analysis_2>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"list"; 
}
