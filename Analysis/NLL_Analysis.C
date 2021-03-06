#include "AddOns/Analysis/Analyses/Analysis_Base.H"

#include "Analysis/Observable_Base.H"
#include "Main/Resum.H"
#include "Analysis/ChannelAlgorithms/ChannelAlgorithm_Base.H"

namespace RESUM {

  class NLL_Analysis: public ANALYSIS::Analysis_Base {
  private:

    typedef std::vector<std::pair<double,double>> CumDist;
    typedef std::map<std::string,CumDist> ObsDist;
    typedef std::map<std::string,ObsDist> ChDist;
    
    ANALYSIS::Argument_Matrix m_params;

    std::vector<Observable_Base*> m_obss;

    Resum *p_resum;

    size_t m_nObs = 0;
    std::vector<std::string> m_ObsNames;
    double m_Ncount = 0;
    double m_sumW = 0;
    double m_sumW2 = 0;
    ObsDist m_NLL;
    ObsDist m_expLO;
    ObsDist m_expNLO;
    ChDist m_NLL_channels;
    ChDist m_expLO_channels;
    ChDist m_expNLO_channels;
    std::vector<std::string> m_channels;

    void Output(const std::string& name, const CumDist& dist,
                const std::vector<double>& xvals);

    std::vector<ChannelAlgorithm_Base::Ptr> m_channelAlgs;    
  public:

    NLL_Analysis(const ANALYSIS::Argument_Matrix &params);

    void Evaluate(const ATOOLS::Blob_List & blobs,
		  double weight,double ncount) override;
    void Evaluate(double weight, double ncount,int mode) override {}
    void Output(const std::string& pname) override;
    void EndEvaluation(double scale) override;
    Primitive_Observable_Base * Copy() const override;

  };// end of class NLL_Analysis

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

NLL_Analysis::NLL_Analysis(const Argument_Matrix& params):
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
  for(size_t i=1; i<params.size(); i++) {
    std::vector<std::string> ps = {params[i].begin()+1,
                                   params[i].end()};
    if (params[i][0] == "ChAlg") {
      m_channelAlgs.emplace_back(ChAlg_Getter::GetObject(params[i][1],
                                                         {params[i][1],
                                                          ps}));
      continue;
    }
    // to stay consistent with previous defaults
    if(m_channelAlgs.size()==0 &&
     reader.GetValue<std::string>("RESUM::SORT_CHANNELS","YES")!="NO") {
      int nborn = p_resum->nLegs()-2;
      m_channelAlgs.emplace_back(ChAlg_Getter::GetObject("KT2_ee",
                                                         {"KT2_ee",
                                                          {std::to_string(nborn),
                                                           ""}}));
    }

    if (params[i].size()<5) continue;
    m_nObs++;
    const double xmin = ATOOLS::ToType<double>(ip->Interprete(params[i][1]));
    const double xmax = ATOOLS::ToType<double>(ip->Interprete(params[i][2]));
    const int nbin = ATOOLS::ToType<int>(ip->Interprete(params[i][3]));
    const string& htype = params[i][4];
    msg_Debugging()<<"Going to add "<<params[i][0]<<" xmin = "<<xmin
                   <<" xmax = "<<xmax<<" Nbin = "<<nbin<<" Type = "<<htype<<"\n";
    Histogram dummy = {HistogramType(htype),xmin,xmax,nbin,"Dummy"};
    vector<double> xvals(dummy.Nbin()+1);
    xvals[0] = dummy.LowEdge(0);
    for(size_t i=0; i<dummy.Nbin(); i++) xvals[i+1] = dummy.HighEdge(i);

    // returns index of observable
    const string& obsName = p_resum->AddObservable({params[i][0],ps},xvals);
    m_ObsNames.push_back(obsName);
    msg_Debugging()<<"Add "<<params[i][0]<<" as "<<obsName<<"\n";
    m_NLL[obsName] = CumDist(xvals.size(),{0,0});
    m_expLO[obsName] = CumDist(xvals.size(),{0,0});
    m_expNLO[obsName] = CumDist(xvals.size(),{0,0});    
    
    
  }
}

void NLL_Analysis::Evaluate(const ATOOLS::Blob_List& blobs,
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
  for(const string& channel: channels) {
    msg_Debugging()<<"Found channel "<< channel<<"...\n";
    if(std::find(m_channels.begin(), m_channels.end(), channel)
       == m_channels.end()) {
      msg_Debugging()<<"... encountered for first time, initalize.\n";
      m_channels.push_back(channel);
      m_NLL_channels[channel] = ObsDist();
      m_expLO_channels[channel] = ObsDist();
      m_expNLO_channels[channel] = ObsDist();
      for(const std::string& name: m_ObsNames) {
        m_NLL_channels[channel][name] = CumDist(m_NLL[name].size(),{0,0});
        m_expLO_channels[channel][name] = CumDist(m_expLO[name].size(),{0,0});
        m_expNLO_channels[channel][name] = CumDist(m_expNLO[name].size(),{0,0});
      }
    }
    else {
      msg_Debugging()<<" ... recognized from before.\n";
    }
  }
  // get results from resum and fill
  for(size_t k=0; k<m_nObs; k++) {
    const std::string& name = m_ObsNames[k];
    const std::vector<double>& resNLL = p_resum->resNLL(k);
    const std::vector<double>& resExpLO = p_resum->resExpLO(k);
    const std::vector<double>& resExpNLO = p_resum->resExpNLO(k);
    const size_t n = p_resum->m_xvals[k].size();
    for(size_t i=0; i<n; i++) {
      const double nll = weight*resNLL[i];
      const double expLO = weight*resExpLO[i];
      const double expNLO = weight*resExpNLO[i];
      msg_Debugging()<<name<<" "<<p_resum->m_xvals[k][i]<<" "<<nll<<" "<<expLO<<" "<<expNLO<<".\n";
      m_NLL[name][i].first += nll;
      m_NLL[name][i].second += sqr(nll);
      m_expLO[name][i].first += expLO;
      m_expLO[name][i].second += sqr(expLO);
      m_expNLO[name][i].first += expNLO;
      m_expNLO[name][i].second += sqr(expNLO);

      for(const string& channel: channels) {
        m_NLL_channels[channel][name][i].first += nll;
        m_NLL_channels[channel][name][i].second += sqr(nll);
        m_expLO_channels[channel][name][i].first += expLO;
        m_expLO_channels[channel][name][i].second += sqr(expLO);
        m_expNLO_channels[channel][name][i].first += expNLO;
        m_expNLO_channels[channel][name][i].second += sqr(expNLO);
      }
    }
  }
}


void NLL_Analysis::EndEvaluation(double scale) {
  for(size_t k=0; k<m_nObs; k++) {
    const std::string& name = m_ObsNames[k];
    const size_t n = p_resum->m_xvals[k].size();
    for(size_t i=0; i<n; i++) {
      m_NLL[name][i].first /= m_Ncount;
      m_NLL[name][i].second /= m_Ncount;
      m_expLO[name][i].first /= m_Ncount;
      m_expLO[name][i].second /= m_Ncount;
      m_expNLO[name][i].first /= m_Ncount;
      m_expNLO[name][i].second /= m_Ncount;

      for(const std::string& ch: m_channels) {
        m_NLL_channels[ch][name][i].first /= m_Ncount;
        m_NLL_channels[ch][name][i].second /= m_Ncount;
        m_expLO_channels[ch][name][i].first /= m_Ncount;
        m_expLO_channels[ch][name][i].second /= m_Ncount;
        m_expNLO_channels[ch][name][i].first /= m_Ncount;
        m_expNLO_channels[ch][name][i].second /= m_Ncount;
      }
    }
  }  
}

void NLL_Analysis::Output(const std::string& name, const CumDist& dist,
                          const std::vector<double>& xvals) {
  msg_Debugging()<<"Writing to "<<name<<"\n";
  
  My_Out_File ofile(name);
  ofile.Open();
  ofile->precision(ToType<int>(rpa->gen.Variable("HISTOGRAM_OUTPUT_PRECISION")));
  *ofile<<"# v Sigma{sumW sumW2 err} barSigma{sumW sumW2 err} NumEntries\n";
  const double sigma = dist[dist.size()-1].first;
  const double sigma2 = dist[dist.size()-1].second;
  for(size_t i=0; i<xvals.size(); i++) {
    const double sig = dist[i].first;
    const double sig2 = dist[i].second;
    const double err = sqrt(sig2-sqr(sig))/(m_Ncount-1);
    const double Bsig = sigma - dist[i].first;
    const double Bsig2 = sigma2 - dist[i].second;
    const double Berr = sqrt(sig2-sqr(sig) + sigma2-sqr(sigma))/(m_Ncount-1);
    *ofile<<xvals[i]<<" "<<sig<<" "<<sig2<<" "<<err<<" ";
    *ofile<<Bsig<<" "<<Bsig2<<" "<<Berr<<" ";
    *ofile<<m_Ncount<<"\n";
  }
}
  
void NLL_Analysis::Output(const std::string& pname) {
  MakeDir(pname+"/NLL");
  MakeDir(pname+"/LO");
  MakeDir(pname+"/NLO");
  for(size_t k=0; k<m_nObs; k++) {
    const std::string& name = m_ObsNames[k];
    Output(pname+"/NLL/"+name+"_Sigma.dat",m_NLL[name],p_resum->m_xvals[k]);
    Output(pname+"/LO/"+name+"_Sigma.dat",m_expLO[name],p_resum->m_xvals[k]);
    Output(pname+"/NLO/"+name+"_Sigma.dat",m_expNLO[name],p_resum->m_xvals[k]);
    for(const std::string& ch: m_channels) {
      Output(pname+"/NLL/Channel_"+ch+"_"+name+"_Sigma.dat",
             m_NLL_channels[ch][name], p_resum->m_xvals[k]);
      Output(pname+"/LO/Channel_"+ch+"_"+name+"_Sigma.dat",
             m_expLO_channels[ch][name], p_resum->m_xvals[k]);
      Output(pname+"/NLO/Channel_"+ch+"_"+name+"_Sigma.dat",
             m_expNLO_channels[ch][name], p_resum->m_xvals[k]);
    }
  }  
}


Primitive_Observable_Base *NLL_Analysis::Copy() const 
{
  return new NLL_Analysis(m_params);
}

DECLARE_GETTER(NLL_Analysis,"NLLResum",Primitive_Observable_Base,Argument_Matrix);

Primitive_Observable_Base *ATOOLS::Getter
<Primitive_Observable_Base,Argument_Matrix,NLL_Analysis>::
operator()(const Argument_Matrix &parameters) const
{
  if (parameters.size()==0 || parameters[0].size()<1) return nullptr;
  return new NLL_Analysis(parameters);
}

void ATOOLS::Getter
<Primitive_Observable_Base,Argument_Matrix,NLL_Analysis>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"list"; 
}
