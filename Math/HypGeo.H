#ifndef _HYPGEO_H_INCLUDED_
#define _HYPGEO_H_INCLUDED_


double _HypGeo_3F2_Py(double a1, double a2, double a3, double b1, double b2, double x);
double _HypGeo_3F2_Arb(double a1, double a2, double a3, double b1, double b2, double x,
                       int regularized=0, int prec=64);

double HypGeo_3F2(double a1, double a2, double a3, double b1, double b2, double x) {
#ifdef USING_ARBLIB
  return _HypGeo_3F2_Arb(a1,  a2,  a3,  b1,  b2,  x);
#elif USING_PYTHON
  return _HypGeo_3F2_Py(a1,  a2,  a3,  b1,  b2,  x);
#else
  THROW(fatal_error,"No hypergeometric function available. Use --enable-arblib or --enable-python.");
  return 0;
#endif
}

#endif
