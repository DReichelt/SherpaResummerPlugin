#ifndef RESUM_THRUSTMAJOR_H
#define RESUM_THRUSTMAJOR_H

#include "Analysis/Observable_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "Observables/Algorithms/ThrustFinder.H"
#include "FFunction/FFunctions.H"

namespace RESUM {
  template <typename T>
  class ThrustMajor_Template: public Observable_Base_Template<T> {
  private:
    using Vec4T = ATOOLS::Vec4<T>;
    using Vec3T = ATOOLS::Vec3<T>;
    
  public:
    ThrustMajor_Template(const Observable_Key& args) : Observable_Base_Template<T>(args) { }

    Obs_Params Parameters(const std::vector<Vec4T>& moms,
                          const std::vector<ATOOLS::Flavour>& flavs,
                          const size_t& l=0) {
      return {1,0,log(2),0};
    }

    using Observable_Base_Template<T>::Name;
    
    T Value(const std::vector<Vec4T>& moms,
            const std::vector<ATOOLS::Flavour>& flavs,
            const size_t& nin) {
      std::map<std::string, typename Algorithm<T>::Ptr> dummy;
      return Value(moms, flavs, dummy, nin);
    }

    std::function<double(double,double&)> FFunction(const std::vector<ATOOLS::Vec4D>& p,
                                                    const std::vector<ATOOLS::Flavour>& fl, 
                                                    const RESUM::Params& params) {
      if(p_F == nullptr) {
        p_F.reset(new FFUNCTION::FFunction(Name()+".dat",-0.468));
      }
      return *p_F;
    }
       
    T Value(const std::vector<Vec4T>& moms,
            const std::vector<ATOOLS::Flavour>& flavs,
            std::map<std::string, typename Algorithm<T>::Ptr>& algorithms,
            const size_t& nin) {
      const std::string name = "ThrustFinder";
      auto alg = algorithms.find(name);
      if(alg==algorithms.end()) {
        alg = algorithms.insert({name,
                                 GetAlgorithm<T>(name, moms, flavs, nin)}).first;
        msg_Debugging()<<"Found jets.\n";
      }
      else msg_Debugging()<<"Reusing jets found earlier.\n";
      return std::dynamic_pointer_cast<ThrustFinder<T>>(alg->second)->ThrustMajor();
    }

    FFUNCTION::FFunction::Ptr p_F = nullptr;
  };

}


#endif
