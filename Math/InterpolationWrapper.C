#include "ATOOLS/Org/Exception.H"
#include "Math/InterpolationWrapper.H"
#include "Math/spline.hpp"


InterpolationWrapper::InterpolationWrapper(const std::vector<double>& xvals,
                                           std::vector<double> yvals) :
    m_xvals(xvals), m_yvals(yvals), m_coeff(xvals.size()) {
  _calc();
}

void InterpolationWrapper::_calc() {
  m_coeff.resize(m_xvals.size());
  spline_pchip_set(m_xvals.size(),&m_xvals[0],&m_yvals[0],&m_coeff[0]);
}

double InterpolationWrapper::Interpolate(double x)  {
  if(m_mode != MODE::HERMITE) THROW(fatal_error, "Only hermite interpolation implemented");
  double ret;
  spline_pchip_val(m_xvals.size(),&m_xvals[0],&m_yvals[0],&m_coeff[0],1,&x, &ret);
  return ret;
}
