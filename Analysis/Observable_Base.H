#ifndef RESUM__Tools__Observable_Base_H
#define RESUM__Tools__Observable_Base_H

#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "Observables/Algorithms/Algorithm.H"
#include "Tools/StringTools.H"
#include "Main/Resum.H"
#include <functional>

namespace RESUM {
  
  struct Observable_Key {
  public:
    std::string m_name;
    Observable_Key(const std::string &name):
      m_name(name) {}
  };

  struct Obs_Params {
    double m_a, m_b, m_logdbar, m_g;
    Obs_Params() = default;
    inline Obs_Params(const double &a,const double &b,
		      const double &logdbar,const double &g):
      m_a(a), m_b(b), m_logdbar(logdbar), m_g(g) {}
  };

  template <typename T>
  class Observable_Base_Template {
  public:
    typedef ATOOLS::Vec4<T> Vec4T;
    typedef ATOOLS::Vec3<T> Vec3T;
  protected:

    std::string m_name;

    double m_logfac = 1.;
    double m_logpow = 1.;
    double m_endpoint = 1.;

    enum EP_MODE {
                  FIXED = 0,
                  DYNAMIC = 1   << 0,
    };

    EP_MODE m_epmode = EP_MODE::FIXED;
    
    const std::map<std::string,EP_MODE> m_EPModeToEnum = {{"FIXED", EP_MODE::FIXED},
                                                          {"DYNAMIC", EP_MODE::DYNAMIC}};
    
  public:
    typedef std::shared_ptr<Observable_Base_Template> Ptr;
    Observable_Base_Template(const Observable_Key &args) : m_name(args.m_name) {
      ATOOLS::Data_Reader read(" ",";","#","=");
      m_logfac = read.GetValue<double>("RESUM::LOGFAC",1.);
      m_logpow = read.GetValue<double>("RESUM::LOGPOW",1.);
      m_endpoint = read.GetValue<double>("RESUM::"+m_name+"::ENDPOINT",1.);
      const std::string& ep_mode = read.GetValue<std::string>("RESUM::"+m_name+"::EP_MODE","FIXED");
      if(is_int(ep_mode)) {
        m_epmode = static_cast<EP_MODE>(to_type<int>(ep_mode));
      }
      else {
        for(const std::string& m: split(ep_mode,"\\|")) {
          m_epmode = static_cast<EP_MODE>(m_epmode | m_EPModeToEnum.at(m));
        }
      }
      
    }

    virtual ~Observable_Base_Template() {};

    // deprecated
    // Obs_Params Parameters
    // (const Vec4T *p,const ATOOLS::Flavour *fl,
    //  const size_t &n,const size_t &l) {
    //   std::vector<Vec4T> moms;
    //   std::vector<ATOOLS::Flavour> flavs;
    //   for(size_t i=0; i<n; i++) {
    //     moms.push_back(p[i]);
    //     flavs.push_back(fl[i]);
    //   }
    //   return Parameters(moms,flavs,l);
    // }

    // // use this
    // virtual Obs_Params Parameters
    // (const std::vector<Vec4T>& p,
    //  const std::vector<ATOOLS::Flavour>& fl,
    //  const size_t &l) = 0;

    Obs_Params Parameters(ATOOLS::Cluster_Amplitude* ampl, 
                          const size_t l) {
      DEBUG_FUNC(this);
      // return this->Parameters(ampl,l);
      std::vector<Vec4T> moms(ampl->Legs().size());
      std::vector<ATOOLS::Flavour> flavs(ampl->Legs().size());
      Vec4T cme = {T(0),T(0),T(0),T(0)};
      for (size_t i=0; i<ampl->Legs().size(); ++i) {
        if(i<ampl->NIn()) cme += -ampl->Leg(i)->Mom();
        moms[i]=i<ampl->NIn()?-ampl->Leg(i)->Mom():ampl->Leg(i)->Mom();
        flavs[i]=i<ampl->NIn()?ampl->Leg(i)->Flav().Bar():ampl->Leg(i)->Flav();
      }
      Obs_Params p = Parameters(moms,flavs,l);
      double logdbar = p.m_logdbar + p.m_a/2.*log(ampl->MuQ2()/cme.Abs2()); 
      return {p.m_a,p.m_b,logdbar,p.m_g};
    }

    virtual Obs_Params Parameters(const std::vector<Vec4T>& p,
                                  const std::vector<ATOOLS::Flavour>& fl,
                                  const size_t& l) = 0;

    virtual std::function<double(double,double&)> FFunction(const std::vector<Vec4T>& p,
                                                            const std::vector<ATOOLS::Flavour>& fl,
                                                            const RESUM::Params& params) = 0;
    
    virtual double Shift(ATOOLS::NLO_subevt *sub) {
        ATOOLS::Vec4D p(sub->p_mom[sub->m_ijt]);
        ATOOLS::Poincare cms(sub->p_mom[0]+sub->p_mom[1]);
        cms.Boost(p);
        std::vector<Vec4T> moms(sub->m_n);
        std::vector<ATOOLS::Flavour> fl(sub->m_n);
        for(int i=0; i<sub->m_n; i++) {
          moms[i] = sub->p_mom[i];
          fl[i] = sub->p_fl[i];
        }
        // Obs_Params ps(Parameters(moms,fl,sub->m_ijt));
        // return ps.m_logdbar-ps.m_b/2.0*log(ATOOLS::sqr(2.0*p[0])/sub->m_mu2[ATOOLS::stp::res]);
        return 0;
    }
    
    double LogArg(double v, ATOOLS::Cluster_Amplitude* ampl) {
      double endpoint = Endpoint(ampl);
      double logfac = LogFac(ampl);
      double logpow = LogPow(ampl);
      return LogArg(v, logfac, logpow, endpoint);
    }

    double ObsVal(double logarg, ATOOLS::Cluster_Amplitude* ampl) {
      double endpoint = Endpoint(ampl);
      double logfac = LogFac(ampl);
      double logpow = LogPow(ampl);
      return ObsVal(logarg, logfac, logpow, endpoint);
    }

    double Endpoint(ATOOLS::Cluster_Amplitude* ampl) {
      if( (m_epmode & EP_MODE::DYNAMIC) != 0) {
        return _Endpoint_Dynamic(ampl);
      }
      return m_endpoint;
    }


    double LogFac() {return m_logfac;}
    double LogPow() {return m_logpow;}

    virtual double LogFac(ATOOLS::Cluster_Amplitude* ampl) {
      double avg = 0;
      double n = 0;
      for(int i=0; i<ampl->Legs().size(); i++) {
        if(ampl->Leg(i)->Flav().StrongCharge() != 0) {
          n++;
          avg += Parameters(ampl,i).m_logdbar;
        }
      }
      avg /= n;
      return exp(avg)*LogFac();
    }
    

    virtual double LogPow(ATOOLS::Cluster_Amplitude* ampl) {
      return LogPow();
    }

    double LogArg(double v, double logfac, 
                  double logpow, double endpoint) {
      double ret = 1;
      if (v < endpoint) {
        ret += pow(logfac/v,logpow) - pow(logfac/endpoint,logpow);
      }
      return 1./pow(ret,1./logpow);
    }

    double ObsVal(double logarg, double logfac, 
                  double logpow, double endpoint) {
      double ret = pow(1./logarg,logpow)+pow(logfac/endpoint,logpow)-1;
      return logfac/pow(ret,1./logpow);
    }

    virtual double _Endpoint_Dynamic(ATOOLS::Cluster_Amplitude* ampl) {
      THROW(not_implemented,"No dynamic endpoint implemented for "+m_name+".");
    }
    
    // deprecated
    T Value
    (const Vec4T *p,const ATOOLS::Flavour *fl,
     const size_t &n,const size_t &nin=2) {
      std::vector<Vec4T> moms;
      std::vector<ATOOLS::Flavour> flavs;
      for(size_t i=0; i<n; i++) {
        moms.push_back(p[i]);
        flavs.push_back(fl[i]);
      }
      return Value(moms, flavs, nin);
    }


    virtual T Value(const std::vector<Vec4T>& p,
                         const std::vector<ATOOLS::Flavour>& fl,
                         const size_t &nin=2) = 0;

    
    T Value
    (const std::vector<Vec4T>& p,
     const std::vector<ATOOLS::Flavour>& fl,
     std::map<std::string, typename Algorithm<T>::Ptr> algorithms,
     const size_t &nin=2) {
      return Value(p, fl, nin);
    }

    
    T Value(ATOOLS::Particle_List *const pl) {
        //TODO: make this work with T
      std::vector<ATOOLS::Vec4D> p(pl->size());
      std::vector<ATOOLS::Flavour> fl(pl->size());
      for (size_t i(0);i<pl->size();++i) {
        p[i]=(*pl)[i]->Momentum();
        fl[i]=(*pl)[i]->Flav();
      }
      return Value(&p.front(),&fl.front(),pl->size(),0);
    }

    inline const std::string Name() { return m_name; }

  };// end of class Observable_Base

  typedef Observable_Base_Template<double> Observable_Base;

  typedef ATOOLS::Getter_Function
  <Observable_Base,Observable_Key> Observable_Getter;
  
}// end of namespace RESUM

#endif


