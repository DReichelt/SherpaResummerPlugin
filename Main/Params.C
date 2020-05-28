#include "Main/Params.H"
#include "MODEL/Main/Running_AlphaS.H"

using namespace RESUM;

double Params::beta0(double scale2) const {
  if(!m_largeNC) return p_as->Beta0(scale2)/M_PI;
  else return 11./12./M_PI;
}
double Params::beta1(double scale2) const {
  if(!m_largeNC) return p_as->Beta1(scale2)/(M_PI*M_PI);
  else return 17./24./M_PI/M_PI;
}
double Params::K_CMW(double scale2) const {
  if(!m_largeNC) return  m_CA*(67./18. - M_PI*M_PI/6.)- m_TR*p_as->Nf(scale2)*10./9.;
  else return 67./18. - M_PI*M_PI/6.;
}

double Params::CollDimGlue(double scale2) const {
  if(!m_largeNC) return -M_PI*beta0(scale2)/m_CA;
  else return -M_PI*beta0(scale2);
}

double Params::CollDimQuark(double scale2) const {
  return -3./4.;
}

double Params::alphaS(double scale2) const {
  return m_constAlpha > 0 ? m_constAlpha : (*p_as)(scale2);
}

