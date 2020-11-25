#ifndef RESUM__READINFUNCTION__READINFUNCTIONS_H
#define RESUM__READINFUNCTION__READINFUNCTIONS_H

#include "ATOOLS/Math/MathTools.H"
#include "Math/InterpolationWrapper.H"
#include <memory>
#include <vector>

namespace RESUM {
    class ReadInFunction : public InterpolationWrapper {
    public:
      typedef std::unique_ptr<ReadInFunction> Ptr;

      ReadInFunction(const std::string& filename, int inc=1,  double yvar=0);
      ReadInFunction(const std::string& filename, double expNLL_LO, 
                     double expNLL_NLO, int inc=1,  double yvar=0);
      
      double operator()(const double x, double& expNLL_NLO);
      double operator()(const double x, double& expNLL_LO, double& expNLL_NLO);

      // void PrintYODA(const std::string& name);

    private:
      void _read(const std::string& filename,int inc, double yvar);
      std::vector<double> m_yerrs;
      double m_expNLL_LO = 0;
      double m_expNLL_NLO = 0;
    };
}

#endif
