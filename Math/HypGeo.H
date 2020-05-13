#ifndef _HYPGEO_H_INCLUDED_
#define _HYPGEO_H_INCLUDED_
/** \file */

/// Hypergeometric function 3F2.
/**
   Function to evaluate the hypergeometric function 3F2.
*/


namespace RESUM {
  /// \name Main Function.
  /// Main function using the optimal implementation available
  double HypGeo_3F2(double a1, double a2, double a3, double b1, double b2, double x);
}

#endif
