#ifndef RESUM__FFUNCTION__FFUNCTIONS_H
#define RESUM__FFUNCTION__FFUNCTIONS_H

namespace RESUM {
  namespace FFUNCTION {
    inline double Additive(const double& Rp) {
      return exp(-GAMMA_E*Rp-ATOOLS::Gammln(1.+Rp));
    }
  }
}

#endif
