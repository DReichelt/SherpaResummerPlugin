#ifndef RESUM_TOTBROAD_H
#define RESUM_TOTBROAD_H

#include "Analysis/Observable_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "Observables/Algorithms/GetAlgorithm.H"
#include "FFunction/FFunctions.H"

namespace RESUM {
  template <typename T>
  class TotBroad_Template: public Observable_Base_Template<T> {
  private:
    using Vec4T = ATOOLS::Vec4<T>;
    using Vec3T = ATOOLS::Vec3<T>;

  public:
    TotBroad_Template(const Observable_Key& args) : Observable_Base_Template<T>(args) {}

    Obs_Params Parameters(const std::vector<Vec4T>& moms,
                          const std::vector<ATOOLS::Flavour>& flavs,
                          const size_t& l=0) {
      return {1.,0.,0.,1.};
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
        p_F.reset(new FFUNCTION::FFunction(Name()+".dat",-0.6514600236711573));
      }
      return *p_F;
    }

    virtual std::set<size_t> ResumMult() {return {2};}
    virtual size_t ResumQCDorderBorn() {return 0;};
    virtual size_t ResumQCDorderLO() {return 1;}
    virtual size_t ResumQCDorderNLO() {return 2;}

       
    T Value(const std::vector<Vec4T>& moms,
            const std::vector<ATOOLS::Flavour>& flavs,
            std::map<std::string, typename Algorithm<T>::Ptr>& algorithms,
            const size_t& nin) {
            // TODO: actually implment anything but the default...
      const std::string& name = ATOOLS::rpa->gen.Variable("RESUM::HeavyJetMass_JetDefinition","ThrustFinder");
      // TODO: this is not really safe in case there are analyses with different template arguments in place
      // This is currently not done, so it is rather a formal problem, but still something to look out for
      auto alg = algorithms.find(name);
      if(alg==algorithms.end()) {
        alg = algorithms.insert({name,GetAlgorithm<T>(name, moms, flavs, nin)}).first;
        msg_Debugging()<<"Found jets.\n";
      }
      else msg_Debugging()<<"Reusing jets found earlier.\n";
      T broad1 = T(0);
      T broad2 = T(0);
      for(const T& broad: alg->second->jetBroads()) {
        if(broad > broad1) {
          broad2 = broad1;
          broad1 = broad;
        }
        else if(broad > broad2) {
          broad2 = broad;
        }
      }
      return broad1+broad2;    
    }

    FFUNCTION::FFunction::Ptr p_F = nullptr;

  };

}


#endif
