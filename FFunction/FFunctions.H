#ifndef RESUM__FFUNCTION__FFUNCTIONS_H
#define RESUM__FFUNCTION__FFUNCTIONS_H

#include "ATOOLS/Math/MathTools.H"
#include "Math/InterpolationWrapper.H"
#include <memory>
#include <vector>

namespace RESUM {
  namespace FFUNCTION {
    inline double Additive(const double& Rp, double& FexpNLL_NLO) {
      FexpNLL_NLO = -pow(M_PI,2)/12; 
      return exp(-GAMMA_E*Rp-ATOOLS::Gammln(1.+Rp));
    }


    class FFunction : public InterpolationWrapper {
    public:
      typedef std::unique_ptr<FFunction> Ptr;

      FFunction(const std::string& filename);
      
      double operator()(const double Rp, double& FexpNLL_NLO);

      // double LaplacePol(const double Rp) const;

    private:
      std::vector<double> m_Rps;
      std::vector<double> m_Fs;
      // std::vector<double> m_bweights;
    };


  }
}

#endif
