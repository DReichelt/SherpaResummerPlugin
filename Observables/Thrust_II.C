#include "Analysis/Observable_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_Limits.H"
#include "FFunction/FFunctions.H"
#include <vector>       
#include <algorithm>       

#define s_ymax std::numeric_limits<double>::max()

using namespace ATOOLS;

namespace RESUM {

  class Thrust_II: public Observable_Base {
  public:

    Thrust_II(const Observable_Key &args): 
    Observable_Base(args) {}

    Obs_Params Parameters
      (const std::vector<Vec4D>& p,
       const std::vector<Flavour>& fl,
       const size_t& l) {
      return Obs_Params(1.0,1.0,0.0,0.0);
    }
    
    std::function<double(double)> FFunction(const std::vector<ATOOLS::Vec4D>& p,
                                            const std::vector<ATOOLS::Flavour>& fl) {
      return FFUNCTION::Additive;
    }

    
    double Value(const std::vector<Vec4D>& p,
                 const std::vector<Flavour>&  fl,
		 const size_t &nin)
    {
      size_t n = p.size();
      Vec4D Q;
      for (size_t i(nin);i<n;++i)
	if (!Flavour(kf_jet).Includes(fl[i])) Q+=p[i];
      double tauB(0.0), eY(sqrt(Q.PPlus()/Q.PMinus())), Y(log(eY));
      for (size_t i(nin);i<n;++i)
	if (Flavour(kf_jet).Includes(fl[i])) {
	  double pp(dabs(p[i].PPlus())), pm(dabs(p[i].PMinus()));
	  double y((pp>0.0&&pm>0.0)?0.5*log(pp/pm):(pp>0.0?s_ymax:-s_ymax));
	  tauB+=(y>Y)?p[i].PMinus()*eY:p[i].PPlus()/eY;
	}
      tauB/=Q.Mass();
      return tauB;
    }

  };// end of class Thrust_II

}// end of namespace RESUM

using namespace RESUM;

DECLARE_GETTER(Thrust_II,"Thrust_II",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,Thrust_II>::
operator()(const Parameter_Type &args) const 
{ return new Thrust_II(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,Thrust_II>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"Thrust_II"; }
