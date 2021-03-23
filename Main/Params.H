#ifndef RESUM_PARAMS_H
#define RESUM_PARAMS_H

#include "ATOOLS/Org/Message.H"

namespace MODEL {
  class Running_AlphaS;
  class One_Running_AlphaS;
}

namespace PDF {
  class PDF_Base;
}


namespace RESUM {
  struct Params {
  private:
    double m_NC = 3;
    double m_CF = (m_NC*m_NC-1.)/(2.*m_NC);
    double m_CA = m_NC;
    double m_TR = 1./2.;
    bool m_largeNC = false;
    double m_constAlpha = -1;
    MODEL::One_Running_AlphaS* p_as;
    PDF::PDF_Base* p_pdf1;
    PDF::PDF_Base* p_pdf2;
    double m_muR2fac = 1;
    double m_muF2fac = 1;
  public:
    Params() = default;
    Params(MODEL::One_Running_AlphaS *as, PDF::PDF_Base *pdf1, PDF::PDF_Base *pdf2, 
           bool largeNC=false, double constAlpha=-1, double muR2fac=1, double muF2fac=1);
    //   : m_largeNC(largeNC), m_constAlpha(constAlpha), m_muR2fac(muR2fac), m_muF2fac(muF2fac); {
    //   p_as = as;
    //   p_pdf1 = pdf1;
    //   p_pdf2 = pdf2;
    //   if(m_largeNC) {
    //     m_CF = m_TR*m_NC;
    //     m_CA = m_NC;
    //     msg_Debugging()<<"Calculate QCD parameters in large NC limit.\n";
    //   }
    //   else msg_Debugging()<<"Calculate QCD parameters with NC = 3.\n";
    // }
    Params(MODEL::Running_AlphaS *as, PDF::PDF_Base *pdf1, PDF::PDF_Base *pdf2, 
           bool largeNC=false, double constAlpha=-1, double muR2fac=1, double muF2fac=1);

    // Params(MODEL::Running_AlphaS *as, PDF::PDF_Base *pdf1, PDF::PDF_Base *pdf2, 
    //        bool largeNC=false, double constAlpha=-1, double muR2fac=1, double muF2fac=1) 
    //   : m_largeNC(largeNC), m_constAlpha(constAlpha), m_muR2fac(muR2fac), m_muF2fac(muF2fac) {
    //   p_as = as->GetAs();
    //   p_pdf1 = pdf1;
    //   p_pdf2 = pdf2;
    //   if(m_largeNC) {
    //     m_CF = m_TR*m_NC;
    //     m_CA = m_NC;
    //     msg_Debugging()<<"Calculate QCD parameters in large NC limit.\n";
    //   }
    //   else msg_Debugging()<<"Calculate QCD parameters with NC = 3.\n";
    // }

    double NC() const {return m_NC;}
    double CA() const {return m_CA;}
    double CF() const {return m_CF;}
    double TR() const {return m_TR;}
    bool largeNC() const {return m_largeNC;}
    double constAlpha() const {return m_constAlpha;}
    double MuR2fac() const {return m_muR2fac;}
    double MuF2fac() const {return m_muF2fac;}
    double beta0(double scale2) const;
    double beta1(double scale2) const;
    double K_CMW(double scale2) const;
    double CollDimGlue(double scale2) const;
    double CollDimQuark(double scale2) const;
    double alphaS(double scale2) const;
    PDF::PDF_Base* PDF(size_t i);
  };
}

#endif
