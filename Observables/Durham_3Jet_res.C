#include "Analysis/Observable_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/Poincare.H"
#include <vector>       
#include <algorithm>       
#include <assert.h> 

using namespace ATOOLS;

namespace RESUM {

  class Durham_3Jet_res: public Observable_Base {
  public:

    Durham_3Jet_res(const Observable_Key &args): 
    Observable_Base(args) {}

    Obs_Params Parameters
      (const std::vector<ATOOLS::Vec4D>& p,
       const std::vector<ATOOLS::Flavour>& fl,
       const size_t &l) {
      return Obs_Params(2.0,0.0,0.0,0.0);
    }

    std::function<double(double)> FFunction(const std::vector<ATOOLS::Vec4D>& p,
                                            const std::vector<ATOOLS::Flavour>& fl) {
      return [](double Rp) { return exp(-GAMMA_E*Rp-2.58*Gammln(1.+0.385*Rp)); };
    }
    

    double Value(const std::vector<Vec4D>& ip,
                 const std::vector<Flavour>& fl,
		 const size_t &nin)
    { return 1;
    }

  };

}

using namespace RESUM;

DECLARE_GETTER(Durham_3Jet_res,"Durham_3Jet_res",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,Durham_3Jet_res>::
operator()(const Parameter_Type &args) const 
{ return new Durham_3Jet_res(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,Durham_3Jet_res>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"Durham_3Jet_res"; }
