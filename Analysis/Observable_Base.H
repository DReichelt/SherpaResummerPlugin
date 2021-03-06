#ifndef RESUM__Tools__Observable_Base_H
#define RESUM__Tools__Observable_Base_H

#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Math/Poincare.H"
#include "Observables/Algorithms/Algorithm.H"
#include "Tools/StringTools.H"
#include "Tools/Key_Base.H"

#include <limits>
#include <functional>

namespace RESUM {

  struct Params;
  
  struct Observable_Key : public Key_Base {
    using Key_Base::Key_Base;
  };

  struct Obs_Params {
    double m_a;
    double m_b;
    double m_logdbar;
    double m_g;
    int m_deltad;
    double m_etamin = 0;
    Obs_Params() = default;
    inline Obs_Params(const double &a,const double &b,
                      const double &logdbar,const double &g, const double &etamin=0., int deltad = 1):
      m_a(a), m_b(b), m_logdbar(logdbar), m_g(g), m_etamin(etamin), m_deltad(deltad) { }
  };

  enum GROOM_MODE {
                   NONE    = 0,       // no grooming
                   SD      = 1 << 0,  // do sd grooming
                   SD_SOFT = 1 << 1,  // expect non-trivial groomed soft
                   SD_PDF  = 1 << 2,  // expect non-trivial groomed pdf
                   SD_COLL = 1 << 3   // expect non-trivial groomed coll 
  };

  
  template <typename T>
  class Observable_Base_Template {
  public:
    typedef ATOOLS::Vec4<T> Vec4T;
    typedef ATOOLS::Vec3<T> Vec3T;
  protected:

    std::string m_name;
    std::string m_tag;
    
    double m_logfac = 1.;
    double m_logpow = 1.;
    double m_endpoint = 1.;
    bool m_unitNorm = true;

    enum EP_MODE {
                  FIXED = 0,
                  DYNAMIC = 1   << 0,
    };

    EP_MODE m_epmode = EP_MODE::FIXED;
    
    const std::map<std::string,EP_MODE> m_EPModeToEnum = {{"FIXED", EP_MODE::FIXED},
                                                          {"DYNAMIC", EP_MODE::DYNAMIC}};

    const std::map<std::string,GROOM_MODE> m_GroomModeToEnum = {{"NONE", GROOM_MODE::NONE},
                                                                {"SD", GROOM_MODE::SD}};
    GROOM_MODE m_gmode = GROOM_MODE::NONE;
    double m_zcut = 0;
    double m_beta = 0;
    double m_R0 = 1.;
    
  public:
    typedef std::shared_ptr<Observable_Base_Template> Ptr;
    Observable_Base_Template(const Observable_Key &args) : m_name(args.m_name) {
      m_tag = args.KwArg("tag",m_name);
      DEBUG_FUNC(m_tag);
      ATOOLS::Data_Reader read(" ",";","#","=");
      ATOOLS::Algebra_Interpreter* ip = read.Interpreter();
      m_logfac = read.GetValue<double>("RESUM::LOGFAC",
                 to_type<double>(ip->Interprete(args.KwArg("logfac","1"))));
      m_logpow = read.GetValue<double>("RESUM::LOGPOW",
                 to_type<double>(ip->Interprete(args.KwArg("logpow","1"))));
      m_unitNorm = bool(read.GetValue<int>("RESUM::UNITNORM",
                        to_type<int>(args.KwArg("unitnorm","1"))));
      m_endpoint = read.GetValue<double>("RESUM::"+m_name+"::ENDPOINT",
                   to_type<double>(ip->Interprete(args.KwArg("ep","1"))));
      const std::string& ep_mode = read.GetValue<std::string>("RESUM::"+m_name+"::EP_MODE",
                                                              args.KwArg("ep_mode","FIXED"));
      if(is_int(ep_mode)) {
        m_epmode = static_cast<EP_MODE>(to_type<int>(ep_mode));
      }
      else {
        for(const std::string& m: split(ep_mode,"\\|")) {
          m_epmode = static_cast<EP_MODE>(m_epmode | m_EPModeToEnum.at(m));
        }
      }
      /*
      const std::string& gmode = read.GetValue<std::string>("RESUM::"+m_name+"::GROOM","NONE");
      if(is_int(gmode)) {
        m_gmode = static_cast<GROOM_MODE>(to_type<int>(gmode));
      }
      else {
        if(!(gmode=="DEFAULT")) {
          for(const std::string& m: split(gmode,"\\|")) {
            m_gmode = static_cast<GROOM_MODE>(m_gmode | m_GroomModeToEnum.at(m));
          }
        }
      }*/
      /*
      if(m_gmode & GROOM_MODE::SD) {
        m_zcut = read.GetValue<double>("RESUM::"+m_name+"::GROOM::ZCUT",
                                       read.GetValue<double>("RESUM::GROOM::ZCUT", 0.1) );
        m_beta = read.GetValue<double>("RESUM::"+m_name+"::GROOM::BETA",
                                       read.GetValue<double>("RESUM::GROOM::BETA", 0));
        m_R0   = read.GetValue<double>("RESUM::"+m_name+"::GROOM::RAD0",
                                       read.GetValue<double>("RESUM::GROOM::RAD0", 1));
      }*/
    }

    virtual ~Observable_Base_Template() {};


    const std::string& Tag() {return m_tag;}
    

    Obs_Params Parameters(ATOOLS::Cluster_Amplitude* ampl, 
                          const size_t l) {
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
      DEBUG_FUNC(Tag()+" -> "+Name());
      msg_Debugging()<<"Parameters for leg "<<l<<" of "; 
      msg_Debugging()<<*ampl<<"\n\n";
      msg_Debugging()<<"a = "<<p.m_a<<"\n";
      msg_Debugging()<<"b = "<<p.m_b<<"\n";
      msg_Debugging()<<"log(\\bar{d}) = "<<p.m_logdbar<<"\n";
      msg_Debugging()<<"g = "<<p.m_g<<"\n";
      msg_Debugging()<<"\\eta_{min} = "<<p.m_etamin<<".\n";
      return {p.m_a,p.m_b,logdbar,p.m_g,p.m_etamin,p.m_deltad};
    }

    virtual Obs_Params Parameters(const std::vector<Vec4T>& p,
                                  const std::vector<ATOOLS::Flavour>& fl,
                                  const size_t& l) = 0;

    virtual std::function<double(double,double&)> FFunction(const std::vector<Vec4T>& p,
                                                            const std::vector<ATOOLS::Flavour>& fl,
                                                            const RESUM::Params& params) = 0;

    virtual std::function<double(double,double&)> SnglFunction(const std::vector<Vec4T>& p,
                                                               const std::vector<ATOOLS::Flavour>& fl,
                                                               const RESUM::Params& params) {
      return [](double t,double& SnglExpNLL_NLO) {SnglExpNLL_NLO=0.; return 1.; };
    }


    
    virtual double SoftGlobal(ATOOLS::Cluster_Amplitude* ampl, size_t i, size_t j, double scale) {
      if(ampl->Leg(i)->Flav().Strong() and ampl->Leg(j)->Flav().Strong()) {
        ATOOLS::Vec4D pi = (i<ampl->NIn() ? -1:1)*ampl->Leg(i)->Mom();
        ATOOLS::Vec4D pj = (j<ampl->NIn() ? -1:1)*ampl->Leg(j)->Mom();
        return log((pi+pj).Abs()/scale);
      }
      else {
        return 0;
      }
    }

    virtual double SoftNonGlobal(ATOOLS::Cluster_Amplitude* ampl, size_t i, size_t j, double scale) {
      return 0.;
    }

    
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
        Obs_Params ps(Parameters(moms,fl,sub->m_ijt));
        return ps.m_logdbar-ps.m_b/2.0*log(ATOOLS::sqr(2.0*p[0])/sub->m_mu2[ATOOLS::stp::res]);
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
      if(!m_unitNorm) return LogFac();
      double avg = 0;
      double n = 0;
      for(int i=0; i<ampl->Legs().size(); i++) {
        Obs_Params ps = Parameters(ampl,i);
        if(ampl->Leg(i)->Flav().StrongCharge() != 0 and ps.m_deltad > 0) {
          n++;
          avg += ps.m_logdbar;
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


    // stuff for grooming
    virtual GROOM_MODE GroomMode() {
      return m_gmode;
    }
    
    virtual GROOM_MODE GroomMode(double v, ATOOLS::Cluster_Amplitude* ampl,
                                 const size_t &l) {
        if(v > GroomTransitionPoint(ampl, l)) return GROOM_MODE::NONE;
        else return m_gmode;
    }

    T GroomZcut() {
      return m_zcut;
    }

    T GroomBeta() {
      return m_beta;
    }

    virtual double GroomTransitionPoint(ATOOLS::Cluster_Amplitude* ampl,
                                        const size_t &l) {
        if(m_gmode & GROOM_MODE::NONE) return 0.;
        Obs_Params para = Parameters(ampl,l);
        if(para.m_deltad == 0) return 0.;
        std::vector<ATOOLS::Vec4D> p(ampl->Legs().size());
        for (size_t i=0; i<ampl->Legs().size(); ++i) {
            p[i]=i<ampl->NIn()?-ampl->Leg(i)->Mom():ampl->Leg(i)->Mom();
        }
        const ATOOLS::Poincare cms= {p[0]+p[1]};
        ATOOLS::Vec4D pl = p[l];
        ATOOLS::Poincare(p[0]+p[1]).Boost(pl);
        
        double El = pl[0];
        double Q12 = (p[0]+p[1]).Abs();
        
        double logfac = LogFac(ampl);
        
        double etamin = log(2.*El/Q12);
        if(para.m_etamin>0.) {
            etamin+=para.m_etamin-log(2.*El/Q12);
        }
        double Q = sqrt(ampl->MuQ2());
        
        double zcp = pow(2.*cosh(pl.Eta())/m_R0,m_beta)*m_zcut;
        
        return pow(zcp*2.*El/Q, para.m_a)*exp(-etamin*(para.m_b+(1.+m_beta)*para.m_a))*exp(para.m_logdbar)/logfac;
    }
    
    virtual double GroomTransitionPoint(ATOOLS::Cluster_Amplitude* ampl){
        double transp = 0.;
        double n = 0.;
        for(size_t i=0; i<ampl->Legs().size(); ++i){
            double transp_i = GroomTransitionPoint(ampl, i);
            if(transp_i>0.){
                n++;
                transp += transp_i;
            }
        }
        if(n>0) return transp/n;
        else return 0.;
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

    
    virtual T Value
    (const std::vector<Vec4T>& p,
     const std::vector<ATOOLS::Flavour>& fl,
     std::map<std::string, typename Algorithm<T>::Ptr>& algorithms,
     const size_t &nin=2) {
      return Value(p, fl, nin);
    }

    virtual bool VetoEvent(const std::vector<Vec4T>& p,
                           const std::vector<ATOOLS::Flavour>& fl,
                           std::map<std::string, typename Algorithm<T>::Ptr>& algorithms,
                           const size_t &nin=2) {
      return false;
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
