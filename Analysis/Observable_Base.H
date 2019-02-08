#ifndef RESUM__Tools__Observable_Base_H
#define RESUM__Tools__Observable_Base_H

#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "Observables/Algorithms/Algorithm.H"


using ATOOLS::sqr;
using ATOOLS::stp;

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

    typedef ATOOLS::Vec4<T> Vec4T;
    typedef ATOOLS::Vec3<T> Vec3T;
  protected:

    std::string m_name;

  public:
    typedef std::shared_ptr<Observable_Base_Template> Ptr;
    Observable_Base_Template(const Observable_Key &args) : m_name(args.m_name) {}

    virtual ~Observable_Base_Template() {};

    // F-function for additiv Obs
    virtual double CalcF(const double Rp){
      return exp(-GAMMA_E*Rp-ATOOLS::Gammln(1.+Rp)); };

    // deprecated
    Obs_Params Parameters
    (const Vec4T *p,const ATOOLS::Flavour *fl,
     const size_t &n,const size_t &l) {
      std::vector<Vec4T> moms;
      std::vector<ATOOLS::Flavour> flavs;
      for(size_t i=0; i<n; i++) {
        moms.push_back(p[i]);
        flavs.push_back(fl[i]);
      }
      return Parameters(moms,flavs,l);};

    // use this
    virtual Obs_Params Parameters
    (const std::vector<Vec4T>& p,
     const std::vector<ATOOLS::Flavour>& fl,
     const size_t &l) = 0;
   
    
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
        return ps.m_logdbar-ps.m_b/2.0*log(sqr(2.0*p[0])/sub->m_mu2[stp::res]);
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
     std::map<std::string, Algorithm::Ptr> algorithms,
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


