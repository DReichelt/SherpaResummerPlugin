#include "AddOns/Analysis/Analyses/Analysis_Base.H"

#include "Analysis/Observable_Base.H"
#include "Main/Resum.H"

namespace RESUM {

  class NLL_Analysis: public ANALYSIS::Analysis_Base {
  private:

    typedef std::vector<std::pair<double,double>> CumDist;
    typedef std::map<std::string,CumDist> ObsDist;
    typedef std::map<std::string,ObsDist> ChDist;
    
    ANALYSIS::Argument_Matrix m_params;

    std::vector<Observable_Base*> m_obss;

    Resum *p_resum;

    size_t nObs = 0;
    std::vector<std::string> ObsNames;
    double N = 0;
    double sumW = 0;
    double sumW2 = 0;
    ObsDist _h_NLL;
    ObsDist _h_expLO;
    ObsDist _h_expNLO;
    ChDist _h_NLL_channels;
    ChDist _h_expLO_channels;
    ChDist _h_expNLO_channels;
    std::vector<std::string> m_channels;

    void Output(const std::string& name, const CumDist& dist,
                const std::vector<double>& xvals);

    
  public:

    NLL_Analysis(const ANALYSIS::Argument_Matrix &params);

    void Evaluate(const ATOOLS::Blob_List & blobs,
		  double weight,double ncount) override;
    void Evaluate(double weight, double ncount,int mode) {}
    void Output(const std::string& pname) override;
    void EndEvaluation(double scale) override;
    Primitive_Observable_Base * Copy() const;

  };// end of class NLL_Analysis

}// end of namespace RESUM

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
  for(size_t i=1; i<params.size(); i++) {
    if (params[i].size()<5) continue;
    nObs++;
    const string& obsName = params[i][0];
    ObsNames.push_back(obsName);
    const double xmin = ATOOLS::ToType<double>(ip->Interprete(params[i][1]));
    const double xmax = ATOOLS::ToType<double>(ip->Interprete(params[i][2]));
    const int nbin = ATOOLS::ToType<int>(ip->Interprete(params[i][3]));
    const string& htype = params[i][4];
    Histogram dummy = {HistogramType(htype),xmin,xmax,nbin,obsName};
    vector<double> xvals(dummy.Nbin()+1);
    xvals[0] = dummy.LowEdge(0);
    for(size_t i=0; i<dummy.Nbin(); i++) xvals[i] = dummy.HighEdge(i);
    _h_NLL[obsName] = CumDist(xvals.size(),{0,0});
    _h_expLO[obsName] = CumDist(xvals.size(),{0,0});
    _h_expNLO[obsName] = CumDist(xvals.size(),{0,0});    
    p_resum->AddObservable(obsName,xvals);
  }
}

void NLL_Analysis::Evaluate(const ATOOLS::Blob_List& blobs,
			double weight, double ncount)
{
  DEBUG_FUNC("");
  if(!p_resum) THROW(fatal_error, "Plugin not loaded.");
  if(!p_resum->Initialized()) return;
  N += ncount;
  sumW += weight;
  sumW2 =+ sqr(weight);
  Particle_List* all = p_ana->GetParticleList(m_listname);
  if (all == nullptr) AddZero(ncount,0);

  // find out channel
  string channel = "";
  for(Particle* p: *all) {
    if(p->Flav().IsGluon()) channel = "g"+channel;
    else channel += "q";
  }

  if(std::find(m_channels.begin(), m_channels.end(), channel) == m_channels.end()) {
    // channel encountered for first time, initalize
    m_channels.push_back(channel);
    _h_NLL_channels[channel] = ObsDist();
    _h_expLO_channels[channel] = ObsDist();
    _h_expNLO_channels[channel] = ObsDist();
    for(std::string& name: ObsNames) {
      _h_NLL_channels[channel][name] = CumDist(_h_NLL[name].size(),{0,0});
      _h_expLO_channels[channel][name] = CumDist(_h_expLO[name].size(),{0,0});
      _h_expNLO_channels[channel][name] = CumDist(_h_expNLO[name].size(),{0,0});
    }
  }
  
  // get results from resum an fill
  for(size_t k=0; k<nObs; k++) {
    const std::string& name = ObsNames[k];
    const std::vector<double>& resNLL = p_resum->m_resNLL[k];
    const std::vector<double>& resExpLO = p_resum->m_resExpLO[k];
    const std::vector<double>& resExpNLO = p_resum->m_resExpNLO[k];
    const size_t n = p_resum->m_xvals[k].size();
    for(size_t i=0; i<n; i++) {
      const double nll = weight*resNLL[i];
      const double expLO = weight*resExpLO[i];
      const double expNLO = weight*resExpNLO[i];
      _h_NLL[name][i].first += nll;
      _h_NLL[name][i].second += sqr(nll);
      _h_NLL_channels[channel][name][i].first += nll;
      _h_NLL_channels[channel][name][i].second += sqr(nll);
      _h_expLO[name][i].first += expLO;
      _h_expLO[name][i].second += sqr(expLO);
      _h_expLO_channels[channel][name][i].first += expLO;
      _h_expLO_channels[channel][name][i].second += sqr(expLO);
      _h_expNLO[name][i].first += expNLO;
      _h_expNLO[name][i].second += sqr(expNLO);
      _h_expNLO_channels[channel][name][i].first += expNLO;
      _h_expNLO_channels[channel][name][i].second += sqr(expNLO);
    }
  }
}


void NLL_Analysis::EndEvaluation(double scale) {
  for(size_t k=0; k<nObs; k++) {
    const std::string& name = ObsNames[k];
    const size_t n = p_resum->m_xvals[k].size();
    for(size_t i=0; i<n; i++) {
      _h_NLL[name][i].first /= N;
      _h_NLL[name][i].second /= N;
      _h_expLO[name][i].first /= N;
      _h_expLO[name][i].second /= N;
      _h_expNLO[name][i].first /= N;
      _h_expNLO[name][i].second /= N;

      for(const std::string& ch: m_channels) {
        _h_NLL_channels[ch][name][i].first /= N;
        _h_NLL_channels[ch][name][i].second /= N;
        _h_expLO_channels[ch][name][i].first /= N;
        _h_expLO_channels[ch][name][i].second /= N;
        _h_expNLO_channels[ch][name][i].first /= N;
        _h_expNLO_channels[ch][name][i].second /= N;
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
  const double sigma = dist[dist.size()-1].first/N;
  const double sigma2 = dist[dist.size()-1].second/N;
  for(size_t i=0; i<xvals.size(); i++) {
    const double sig = dist[i].first;
    const double sig2 = dist[i].second;
    const double err = sqrt(sig2-sqr(sig))/(N-1);
    const double Bsig = sigma - dist[i].first;
    const double Bsig2 = sigma2 - dist[i].second;
    const double Berr = sqrt(sig2-sqr(sig) + sigma2-sqr(sigma))/(N-1);
    *ofile<<xvals[i]<<" "<<sig<<" "<<sig2<<" "<<err<<" ";
    *ofile<<Bsig<<" "<<Bsig2<<" "<<Berr<<"\n";
  }
}
  
void NLL_Analysis::Output(const std::string& pname) {
  MakeDir(pname+"/NLL");
  MakeDir(pname+"/LO");
  MakeDir(pname+"/NLO");
  for(size_t k=0; k<nObs; k++) {
    const std::string& name = ObsNames[k];
    Output(pname+"/NLL/"+name+"_Sigma.dat",_h_NLL[name],p_resum->m_xvals[k]);
    Output(pname+"/LO/"+name+"_Sigma.dat",_h_expLO[name],p_resum->m_xvals[k]);
    Output(pname+"/NLO/"+name+"_Sigma.dat",_h_expNLO[name],p_resum->m_xvals[k]);
    for(const std::string& ch: m_channels) {
      Output(pname+"/NLL/Channel_"+ch+"_"+name+"_Sigma.dat",
             _h_NLL_channels[ch][name], p_resum->m_xvals[k]);
      Output(pname+"/LO/Channel_"+ch+"_"+name+"_Sigma.dat",
             _h_expLO_channels[ch][name], p_resum->m_xvals[k]);
      Output(pname+"/NLO/Channel_"+ch+"_"+name+"_Sigma.dat",
             _h_expNLO_channels[ch][name], p_resum->m_xvals[k]);
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
  if (parameters.size()==0 || parameters[0].size()<1) return NULL;
  return new NLL_Analysis(parameters);
}

void ATOOLS::Getter
<Primitive_Observable_Base,Argument_Matrix,NLL_Analysis>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"list"; 
}
