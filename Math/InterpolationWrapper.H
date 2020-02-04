#ifndef RESUM_INTERPOLATION_WRAPPER_H
#define RESUM_INTERPOLATION_WRAPPER_H
#include <vector>

/// Interpolation wrapper class.
/**
   Wrapper for (hermite) interpolation from linpack.
*/


class InterpolationWrapper {
public:
  /// Default constructor.
  InterpolationWrapper() = default;

  /// Constructor 
  /**
     \param xvals Array of x-values
     \param yvals Array of y-values interpreted of \p yvals[i] = f(\p xvals[i])
  */
  InterpolationWrapper(const std::vector<double>& xvals,
                       std::vector<double> yvals);

  /// Interpolated value f(\p x).
  double Interpolate(double x);

  enum MODE {
    DEFAULT = 0,
    HERMITE = 1
  };
  
protected:

  /// Helper function doing the work in constructor.
  void _calc();

  /// Mode to do interpolation
  MODE m_mode = MODE::DEFAULT;

  /// interpolation grid
  std::vector<double> m_xvals;

  /// values to interpolate
  std::vector<double> m_yvals;

private:
  /// vector of coefficients (internal)
  std::vector<double> m_coeff;
};

#endif
