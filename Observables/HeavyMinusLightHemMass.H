#ifndef RESUM_HEAVYMINUSLIGHTMASS_H
#define RESUM_HEAVYMINUSLIGHTMASS_H

#include "Analysis/Observable_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "Observables/Algorithms/Algorithm.H"
#include "FFunction/FFunctions.H"

namespace RESUM {
  template <typename T>
  class HeavyMinLightMass_Template: public Observable_Base_Template<T> {
  private:
    using Vec4T = ATOOLS::Vec4<T>;
    using Vec3T = ATOOLS::Vec3<T>;

  public:
    HeavyMinLightMass_Template(const Observable_Key& args) : Observable_Base_Template<T>(args) {}

    Obs_Params Parameters(const std::vector<Vec4T>& moms,
                          const std::vector<ATOOLS::Flavour>& flavs,
                          const size_t& l=0) {
      return {1.,1.,0.,1.};
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
       
    T Value(const std::vector<Vec4T>& moms,
            const std::vector<ATOOLS::Flavour>& flavs,
            std::map<std::string, typename Algorithm<T>::Ptr>& algorithms,
            const size_t& nin) {
      const std::string& name = ATOOLS::rpa->gen.Variable("RESUM::HeavyJetMass_JetDefinition","ThrustFinder");
      // TODO: this is not really safe in case there are analyses with different template arguments in place
      // This is currently not done, so it is rather a formal problem, but still something to look out for
      auto alg = algorithms.find(name);
      if(alg==algorithms.end()) {
        alg = algorithms.insert({name,GetAlgorithm<T>(name, moms, flavs, nin)}).first;
        msg_Debugging()<<"Found jets.\n";
      }
      else msg_Debugging()<<"Reusing jets found earlier.\n";
      const auto& jet_inds = alg->second->Jets();
      Vec4T sum = {0,0,0,0};
      for(size_t i=nin; i<moms.size(); i++) {
        sum += moms[i];
      }
      msg_Debugging()<<"Calculated event energy.\n";
      T mass2h = 0;
      T mass2l = 0;
      for(const T& m2: alg->second->jetMasses()) {
        if(m2 > mass2h) {
          mass2l = mass2h;
          mass2h = m2;
        }
        else if(m2 > mass2l) {
          mass2l = m2;
        }
      }
      return (mass2h-mass2l)/alg->second->Q2();
    }

  };

}


#endif
