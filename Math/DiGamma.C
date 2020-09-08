#include "Math/DiGamma.H"
#include "ATOOLS/Org/Exception.H"

// using namespace RESUM;

double _DiGamma_Arb(double x, int prec=64);

namespace RESUM {

  double DiGamma(double x) {
#ifdef USING_ARBLIB
    return _DiGamma_Arb( x);
#else
    THROW(fatal_error,"No digamma function available. Use --enable-arblib.");
    return 0;
#endif
  }

}

#ifdef USING_ARBLIB
#include "arb.h"

double _DiGamma_Arb(double z, int prec){
  // init arb variables
  arb_t result; arb_init(result);
  arb_t arg; arb_init(arg); arb_set_d(arg,z);

  // actual calculation
  arb_digamma(result, arg, prec);

  // clean up
  arb_clear(arg);
  const double ret = arf_get_d(arb_midref(result), ARF_RND_NEAR);
  arb_clear(result);

  return ret;
}
#endif
