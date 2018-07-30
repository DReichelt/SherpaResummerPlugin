#ifndef RESUM__Tools__Observable_Base_H
#define RESUM__Tools__Observable_Base_H

#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "Observables/Algorithms/Algorithm.H"

namespace RESUM {
  
  struct Observable_Key {
  public:
    std::string m_name;
    Observable_Key(const std::string &name):
      m_name(name) {}
  };

  struct Obs_Params {
    double m_a, m_b, m_d, m_g;
    inline Obs_Params(const double &a,const double &b,
		      const double &d,const double &g):
      m_a(a), m_b(b), m_d(d), m_g(g) {}
  };

  class Observable_Base {
  protected:

    std::string m_name;

  public:

    Observable_Base(const Observable_Key &args);

    virtual ~Observable_Base();

    Obs_Params Parameters
    (const ATOOLS::Vec4D *p,const ATOOLS::Flavour *fl,
     const size_t &n,const size_t &l) {
      std::vector<ATOOLS::Vec4D> moms;
      std::vector<ATOOLS::Flavour> flavs;
      for(size_t i=0; i<n; i++) {
        moms.push_back(p[i]);
        flavs.push_back(fl[i]);
      }
      return Parameters(moms,flavs,l);};

    virtual Obs_Params Parameters
    (const std::vector<ATOOLS::Vec4D>& p,
     const std::vector<ATOOLS::Flavour>& fl,
     const size_t &l) = 0;
   
    
    virtual double Shift(ATOOLS::NLO_subevt *sub);

    double Value
    (const ATOOLS::Vec4D *p,const ATOOLS::Flavour *fl,
     const size_t &n,const size_t &nin=2) {
      std::vector<ATOOLS::Vec4D> moms;
      std::vector<ATOOLS::Flavour> flavs;
      for(size_t i=0; i<n; i++) {
        moms.push_back(p[i]);
        flavs.push_back(fl[i]);
      }
      return Value(moms, flavs, nin);
    }


    virtual double Value(const std::vector<ATOOLS::Vec4D>& p,
                         const std::vector<ATOOLS::Flavour>& fl,
                         const size_t &nin=2) = 0;

    
    double Value
    (const std::vector<ATOOLS::Vec4D>& p,
     const std::vector<ATOOLS::Flavour>& fl,
     std::map<std::string, Algorithm::Ptr> algorithms,
     const size_t &nin=2) {
      return Value(p, fl, nin);
    }

    
    double Value(ATOOLS::Particle_List *const pl);

    inline const std::string Name() { return m_name; }
        
  };// end of class Observable_Base
  
  typedef ATOOLS::Getter_Function
  <Observable_Base,Observable_Key> Observable_Getter;

}// end of namespace RESUM

#endif


