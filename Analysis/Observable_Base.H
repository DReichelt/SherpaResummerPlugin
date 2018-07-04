#ifndef RESUM__Tools__Observable_Base_H
#define RESUM__Tools__Observable_Base_H

#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Phys/NLO_Subevt.H"

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

    virtual Obs_Params Parameters
    (const ATOOLS::Vec4D *p,const ATOOLS::Flavour *fl,
     const size_t &n,const size_t &l) = 0;

    virtual double Shift(ATOOLS::NLO_subevt *sub);

    virtual double Value
    (const ATOOLS::Vec4D *p,const ATOOLS::Flavour *fl,
     const size_t &n,const size_t &nin=2) = 0;

    double Value(ATOOLS::Particle_List *const pl);

    inline const std::string Name() { return m_name; }
        
  };// end of class Observable_Base
  
  typedef ATOOLS::Getter_Function
  <Observable_Base,Observable_Key> Observable_Getter;

}// end of namespace RESUM

#endif


