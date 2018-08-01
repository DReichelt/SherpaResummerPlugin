#ifndef _MULTIPLEPRECISION_H_
#define _MULTIPLEPRECISION_H_

#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/gmp.hpp>
#include "ATOOLS/Math/Vec4.H"
#include "ATOOLS/Math/Vec3.H"

using namespace ATOOLS;

namespace mp = boost::multiprecision; 

typedef mp::number<mp::mpfr_float_backend<510> >  mpfloat;

DECLARE_PROMOTE(mpfloat,mpfloat,mpfloat);
typedef Vec4<mpfloat> Vec4MP;
typedef Vec3<mpfloat> Vec3MP;


#endif
