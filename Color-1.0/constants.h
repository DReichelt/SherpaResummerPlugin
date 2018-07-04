#if !defined(MYLIB_CONSTANTS_H)
#define MYLIB_CONSTANTS_H 1

namespace constants {
    const double pi    =    3.1415926535;
    const double Ec    =    0.5772156649;
    const double Ct    =    0.9159655941;
    const double Cf    =    4./3.;
    const double Ca    =    3.   ;
    const double Nc    =    3.   ;
    const double nf    =    5.   ;
    const double aM    =    .11  ;
    const double b0    =    (11*Ca - 2*nf)/(12*pi) ;
    const double b1    =    (17*Ca*Ca - 5*Ca*nf - 3*Cf*nf)/(24*pi*pi) ;
    const double Blq   =    -3./4.;
    const double Blg   =    -b0*pi/Ca;
    const double Kp    =    (Ca*(67./18. - pi*pi/6.) - (5./9.)*nf) ;
}
#endif
