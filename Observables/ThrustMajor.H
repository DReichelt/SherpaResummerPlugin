#ifndef RESUM_HEAVYJETMASS_H
#define RESUM_HEAVYJETMASS_H

#include "Analysis/Observable_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "Observables/Algorithms/ThrustFinder.H"
#include "FFunction/FFunctions.H"

using namespace ATOOLS;


namespace RESUM {
  template <typename T>
  class ThrustMajor_Template: public Observable_Base_Template<T> {
  private:
    
  public:
    ThrustMajor_Template(const Observable_Key& args) : Observable_Base_Template<T>(args) { }

    Obs_Params Parameters(const std::vector<Vec4<T>>& moms,
                          const std::vector<Flavour>& flavs,
                          const size_t& l=0) {
      return {1,0,-log(2),2};
    }

    using Observable_Base_Template<T>::Name;
    
    T Value(const std::vector<Vec4<T>>& moms,
                 const std::vector<Flavour>& flavs,
                 const size_t& nin) {
      std::map<std::string, Algorithm::Ptr> dummy;
      return Value(moms, flavs, dummy, nin);
    }

    std::function<double(double)> FFunction(const std::vector<ATOOLS::Vec4D>& p,
                                            const std::vector<ATOOLS::Flavour>& fl) {
      if(p_F == nullptr) {
        p_F.reset(new FFUNCTION::FFunction(Name()+".dat"));
      }
      return *p_F;
    }
       
    T Value(const std::vector<Vec4<T>>& moms,
                 const std::vector<Flavour>& flavs,
                 std::map<std::string, Algorithm::Ptr>& algorithms,
                 const size_t& nin) {
      // TODO: actually implment anything but the default...
      const std::string& name = "ThrustFinder";
      // TODO: this is not really safe in case there are analyses with different template arguments in place
      // This is currently not done, so it is rather a formal problem, but still something to look out for
      auto alg = algorithms.find(name);
      if(alg==algorithms.end()) {
        alg = algorithms.insert({name,GetAlgorithm<T>(name, moms, flavs, nin)}).first;
        msg_Debugging()<<"Found jets.\n";
      }
      else msg_Debugging()<<"Reusing jets found earlier.\n";
      return std::dynamic_pointer_cast<ThrustFinder<T>>(alg->second)->ThrustMajor();
    }

    FFUNCTION::FFunction::Ptr p_F = nullptr;
  };

}


#endif
