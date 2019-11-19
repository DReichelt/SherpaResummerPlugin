#ifndef RESUM_HEAVYJETMASS_H
#define RESUM_HEAVYJETMASS_H

#include "Analysis/Observable_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "Observables/Algorithms/Algorithm.H"

using namespace ATOOLS;


namespace RESUM {
  template <typename T>
  class HeavyJetMass_Template: public Observable_Base_Template<T> {
  private:

  public:
    HeavyJetMass_Template(const Observable_Key& args) : Observable_Base_Template<T>(args) {}

    Obs_Params Parameters(const std::vector<Vec4<T>>& moms,
                          const std::vector<Flavour>& flavs,
                          const size_t& l=0) {
      return {1,1,0,1};
    }

    T Value(const std::vector<Vec4<T>>& moms,
                 const std::vector<Flavour>& flavs,
                 const size_t& nin) {
      std::map<std::string, Algorithm::Ptr> dummy;
      return Value(moms, flavs, dummy, nin);
    }

    std::function<double(double,double&)> FFunction(const std::vector<ATOOLS::Vec4D>& p,
                                                    const std::vector<ATOOLS::Flavour>& fl, 
                                                    const RESUM::Params& params) {
      //@TODO: find expansion
      return [](double Rp,double& Fexp) {Fexp=1.; return exp(-GAMMA_E*Rp-2*Gammln(1.+0.5*Rp)); };
    }
       
    T Value(const std::vector<Vec4<T>>& moms,
                 const std::vector<Flavour>& flavs,
                 std::map<std::string, Algorithm::Ptr>& algorithms,
                 const size_t& nin) {
      // TODO: actually implment anything but the default...
      const std::string& name = rpa->gen.Variable("RESUM::HeavyJetMass_JetDefinition","ThrustFinder");
      // TODO: this is not really safe in case there are analyses with different template arguments in place
      // This is currently not done, so it is rather a formal problem, but still something to look out for
      auto alg = algorithms.find(name);
      if(alg==algorithms.end()) {
        alg = algorithms.insert({name,GetAlgorithm<T>(name, moms, flavs, nin)}).first;
        msg_Debugging()<<"Found jets.\n";
      }
      else msg_Debugging()<<"Reusing jets found earlier.\n";
      const auto& jet_inds = alg->second->Jets();
      Vec4<T> sum = {0,0,0,0};
      for(size_t i=nin; i<moms.size(); i++) {
        sum += moms[i];
      }
      msg_Debugging()<<"Calculated event energy.\n";
      T mass2 = 0;
      for(const auto& inds: jet_inds) {
        Vec4<T> cur = {0,0,0,0};
        msg_Debugging()<<"\n";
        for(int ind: inds) {
          msg_Debugging()<<moms.at(ind)<<"\n";
          cur += moms.at(ind);
        }
        msg_Debugging()<<"\n";
        T m2 = cur.Abs2();
        // TODO: replace this by max function appropriate for T
        mass2 = m2 > mass2 ? m2 : mass2;
      }
      msg_Debugging()<<"Finished.\n";
      return mass2/sum.Abs2();
    }

  };

}


#endif
