#include "Main/Params.H"
#include "MODEL/Main/Running_AlphaS.H"

using namespace RESUM;

Params::Params(MODEL::One_Running_AlphaS *as, PDF::PDF_Base *pdf1, PDF::PDF_Base *pdf2, 
           bool largeNC, double constAlpha, double muR2fac, double muF2fac) 
      : m_largeNC(largeNC), m_constAlpha(constAlpha), m_muR2fac(muR2fac), m_muF2fac(muF2fac) {
  p_as = as;
  p_pdf1 = pdf1;
  p_pdf2 = pdf2;
  if(m_largeNC) {
    m_CF = m_TR*m_NC;
    m_CA = m_NC;
    msg_Debugging()<<"Calculate QCD parameters in large NC limit.\n";
  }
  else msg_Debugging()<<"Calculate QCD parameters with NC = 3.\n";
}


Params::Params(MODEL::Running_AlphaS *as, PDF::PDF_Base *pdf1, PDF::PDF_Base *pdf2, 
       bool largeNC, double constAlpha, double muR2fac, double muF2fac) 
  : Params(as->GetAs(), pdf1, pdf2, largeNC, constAlpha, muR2fac, muF2fac) {}


double Params::beta0(double scale2) const {
  if(!m_largeNC) return p_as->Beta0(m_muR2fac*scale2)/M_PI;
  else return 11./12./M_PI;
}
double Params::beta1(double scale2) const {
  if(!m_largeNC) return p_as->Beta1(m_muR2fac*scale2)/(M_PI*M_PI);
  else return 17./24./M_PI/M_PI;
}
double Params::K_CMW(double scale2) const {
  if(!m_largeNC) return  m_CA*(67./18. - M_PI*M_PI/6.)- m_TR*p_as->Nf(m_muR2fac*scale2)*10./9.;
  else return 67./18. - M_PI*M_PI/6.;
}

double Params::CollDimGlue(double scale2) const {
  if(!m_largeNC) return -M_PI*beta0(m_muR2fac*scale2)/m_CA;
  else return -M_PI*beta0(scale2);
}

double Params::CollDimQuark(double scale2) const {
  return -3./4.;
}

double Params::alphaS(double scale2) const {
  return m_constAlpha > 0 ? m_constAlpha : (*p_as)(m_muR2fac*scale2);
}

PDF::PDF_Base* Params::PDF(size_t i) {
  if(i==0) return p_pdf1;
  if(i==1) return p_pdf2;
  return nullptr;
}
