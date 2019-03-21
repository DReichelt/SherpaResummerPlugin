#ifndef RESUM__FFUNCTION__FFUNCTIONS_H
#define RESUM__FFUNCTION__FFUNCTIONS_H

#include "ATOOLS/Math/MathTools.H"
#include <memory>
#include <vector>

namespace RESUM {
  namespace FFUNCTION {
    inline double Additive(const double& Rp) {
      return exp(-GAMMA_E*Rp-ATOOLS::Gammln(1.+Rp));
    }


    class FFunction {
    public:
      typedef std::unique_ptr<FFunction> Ptr;

      FFunction(const std::string& filename);
      
      double operator()(const double Rp);
      

    private:
      std::vector<double> Rps;
      std::vector<double> Fs;      
    };


  }
}

#endif
