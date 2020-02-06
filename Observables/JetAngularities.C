#include "Analysis/Observable_Base.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_Limits.H"
#include "FFunction/FFunctions.H"
#include <vector>       
#include <algorithm>       

#include "Observables/Algorithms/Algorithm.H"

using namespace ATOOLS;

namespace RESUM {

  class JetAngularities: public Observable_Base {
  public:

    JetAngularities(const Observable_Key &args): 
    Observable_Base(args) {}

    Obs_Params Parameters
      (const std::vector<ATOOLS::Vec4D>& p,
       const std::vector<ATOOLS::Flavour>& fl,
       const size_t& l) {
      // TODO: this is a dummy
      return Obs_Params(2.0,0.0,0.0,0.0);
    }

    std::function<double(double,double&)> FFunction(const std::vector<ATOOLS::Vec4D>& p,
                                                    const std::vector<ATOOLS::Flavour>& fl, 
                                                    const RESUM::Params& params) {
      // TODO: this is a dummy
      return FFUNCTION::Additive;
    }


    double Value(const std::vector<Vec4D>& moms,
            const std::vector<ATOOLS::Flavour>& flavs,
            const size_t& nin) {
      std::map<std::string, typename Algorithm<double>::Ptr> dummy;
      return Value(moms, flavs, dummy, nin);
    }


    double Value(const std::vector<Vec4D>& ip,
            const std::vector<ATOOLS::Flavour>& fl,
            std::map<std::string, typename Algorithm<double>::Ptr>& algorithms,
            const size_t& nin) {

      if(ip.size() <= nin) return 0;
      msg_Debugging()<<"Start jet angularity.\n";
      // some parameters
      // TODO: make these dynamic, avoid magic matching
      const double R_jet = 0.5;          
      const double pt_min_cut = 0.;        
      const double alpha = 2.;
      
      
      const std::string& name = "FJmaxPTjet";
      auto alg = algorithms.find(name);
      if(alg==algorithms.end()) {
        alg = algorithms.insert({name,GetAlgorithm<double>(name, ip, fl, nin)}).first;
        msg_Debugging()<<"Found jets.\n";
      }
      else msg_Debugging()<<"Reusing jets found earlier.\n";

      std::vector<ATOOLS::Vec4D> constits =  alg->second->apply(ip,0);
      ATOOLS::Vec4D axis(0.,alg->second->jetAxes(0));
      
      
      double lambda = 0.;
      
      for (ATOOLS::Vec4D con: constits) {
        
        double z = con.PPerp() / alg->second->Q2();
        
        double Delta = sqrt(sqr(con.Y()-axis.Y()) +
                            sqr(con.DPhi(axis)));
        
        lambda +=  z * pow(Delta / R_jet, alpha);
      }

      return lambda;
    }           
        
    
  };// end of class Y1_II

}// end of namespace RESUM

using namespace RESUM;


DECLARE_GETTER(JetAngularities,"JetAngularities",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,JetAngularities>::
operator()(const Parameter_Type &args) const 
{ return new JetAngularities(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,JetAngularities>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"JetAngularities"; }
