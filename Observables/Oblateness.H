#ifndef RESUM_OBLATENESS_H
#define RESUM_OBLATENESS_H
#include "Analysis/Observable_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "Observables/Algorithms/GetAlgorithm.H"
#include "FFunction/FFunctions.H"

namespace RESUM {
  template <typename T>
  class Oblateness_Template: public Observable_Base_Template<T> {
  private:
    using Vec4T = ATOOLS::Vec4<T>;
    using Vec3T = ATOOLS::Vec3<T>;

  public:
    Oblateness_Template(const Observable_Key& args) : Observable_Base_Template<T>(args) {}

    Obs_Params Parameters(const std::vector<Vec4T>& moms,
                          const std::vector<ATOOLS::Flavour>& flavs,
                          const size_t& l=0) {
      return {1,1,log(2.),1};
    }

    T Value(const std::vector<Vec4T>& moms,
                 const std::vector<ATOOLS::Flavour>& flavs,
                 const size_t& nin) {
      std::map<std::string, typename Algorithm<T>::Ptr> dummy;
      return Value(moms, flavs, dummy, nin);
    }

    std::function<double(double,double&)> FFunction(const std::vector<ATOOLS::Vec4D>& p,
                                                    const std::vector<ATOOLS::Flavour>& fl, 
                                                    const RESUM::Params& params) {
      /// @TODO insert correct function
      return FFUNCTION::Additive;
    }

    virtual std::set<size_t> ResumMult() {return {2};}
    virtual size_t ResumQCDorderBorn() {return 0;};
    virtual size_t ResumQCDorderLO() {return 1;}
    virtual size_t ResumQCDorderNLO() {return 2;}

       
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
      auto thrust = std::dynamic_pointer_cast<ThrustFinder<T>>(alg->second);
      return thrust->ThrustMajor()-thrust->ThrustMinor();
     

    }

  };

}


#endif
