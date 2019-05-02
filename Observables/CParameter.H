#ifndef RESUM_HEAVYJETMASS_H
#define RESUM_HEAVYJETMASS_H

#include "Analysis/Observable_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "Observables/Algorithms/Algorithm.H"
#include "FFunction/FFunctions.H"

using namespace ATOOLS;


namespace RESUM {
  template <typename T>
  class CParameter_Template: public Observable_Base_Template<T> {
  private:

  public:
    CParameter_Template(const Observable_Key& args) : Observable_Base_Template<T>(args) {}

    Obs_Params Parameters(const std::vector<Vec4<T>>& moms,
                          const std::vector<Flavour>& flavs,
                          const size_t& l=0) {
      return {1,1,log(6),1};
    }

    T Value(const std::vector<Vec4<T>>& moms,
                 const std::vector<Flavour>& flavs,
                 const size_t& nin) {
      std::map<std::string, Algorithm::Ptr> dummy;
      return Value(moms, flavs, dummy, nin);
    }

    std::function<double(double)> FFunction(const std::vector<ATOOLS::Vec4D>& p,
                                            const std::vector<ATOOLS::Flavour>& fl) {
      return FFUNCTION::Additive;
    }
       
    T Value(const std::vector<Vec4<T>>& moms,
            const std::vector<Flavour>& flavs,
            std::map<std::string, Algorithm::Ptr>& algorithms,
            const size_t& nin) {
      const std::vector<Vec4<T>>& moms_FF = {moms.begin()+nin,moms.end()};
      Vec4<T> sum = {T(0),T(0),T(0),T(0)};
      for(const auto& p: moms_FF) sum += p;
      T Cpar = T(3);
      for(const auto& p: moms_FF) {
        // Note: default accuracy 1e-12 to strict apparently.
        if(!IsZero(p.Abs2(),1e-6)) THROW(not_implemented, "C parameter is implemented assuming massless particles. p2 = "+ToString(p.Abs2()));
        for(const auto& q: moms_FF) {
          Cpar -= T(3)/T(2) * pow(p*q,2)/(p*sum)/(q*sum);
        }
      }
      return Cpar;
    }

  };

}


#endif
