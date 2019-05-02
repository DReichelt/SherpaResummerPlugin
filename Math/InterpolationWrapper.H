#ifndef RESUM_INTERPOLATION_WRAPPER_H
#define RESUM_INTERPOLATION_WRAPPER_H
#include <vector>


class InterpolationWrapper {
public:
  InterpolationWrapper() = default;
  InterpolationWrapper(const std::vector<double>& xvals,
                       std::vector<double> yvals);

  double Interpolate(double x);

  enum MODE {
    DEFAULT = 0,
    HERMITE = 1
  };
  
protected:

  void _calc();
  MODE m_mode = MODE::DEFAULT;
  std::vector<double> m_xvals;
  std::vector<double> m_yvals;
  std::vector<double> m_coeff;
};

#endif
