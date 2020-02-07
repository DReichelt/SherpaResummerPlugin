#ifndef RESUM_PARAMS_H
#define RESUM_PARAMS_H

#include "ATOOLS/Org/Message.H"

namespace MODEL {
  class Running_AlphaS;
}


namespace RESUM {
  struct Params {
  private:
    double m_NC = 3;
    double m_CF = (m_NC*m_NC-1.)/(2.*m_NC);
    double m_CA = m_NC;
    double m_TR = 1./2.;
    bool m_largeNC = false;
    MODEL::Running_AlphaS *p_as;
  public:
    Params() = default;
    Params(MODEL::Running_AlphaS *as, bool largeNC=false) : m_largeNC(largeNC) {
      p_as = as;
      if(m_largeNC) {
        m_CF = m_TR*m_NC;
        m_CA = m_NC;
        msg_Debugging()<<"Calculate QCD parameters in large NC limit.\n";
      }
      else msg_Debugging()<<"Calculate QCD parameters with NC = 3.\n";
    }
    double NC() const {return m_NC;}
    double CA() const {return m_CA;}
    double CF() const {return m_CF;}
    double TR() const {return m_TR;}
    bool largeNC() const {return m_largeNC;}
    double beta0(double scale2) const;
    double beta1(double scale2) const;
    double K_CMW(double scale2) const;
    double CollDimGlue(double scale2) const;
    double CollDimQuark(double scale2) const;
  };
}

#endif
