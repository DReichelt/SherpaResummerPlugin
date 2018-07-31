#include "Analysis/Observable_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "Observables/Algorithms/Algorithm.H"

using namespace ATOOLS;


namespace RESUM {
  class HeavyJetMass: public Observable_Base {
  private:

  public:
    HeavyJetMass(const Observable_Key& args) : Observable_Base(args) {}

    Obs_Params Parameters(const std::vector<Vec4D>& moms,
                          const std::vector<Flavour>& flavs,
                          const size_t& l=0) {
      return {1,1,1,1};
    }

    double Value(const std::vector<Vec4D>& moms,
                 const std::vector<Flavour>& flavs,
                 const size_t& nin) {
      std::map<std::string, Algorithm::Ptr> dummy;
      return Value(moms, flavs, dummy, nin);
    }

    
    double Value(const std::vector<Vec4D>& moms,
                 const std::vector<Flavour>& flavs,
                 std::map<std::string, Algorithm::Ptr>& algorithms,
                 const size_t& nin) {
      // TODO: actually implment anything but the default...
      const std::string& name = rpa->gen.Variable("RESUM::HeavyJetMass_JetDefinition","ThrustFinder");
      auto alg = algorithms.find(name);
      if(alg==algorithms.end()) {
        alg = algorithms.insert({name,GetAlgorithm(name, moms, flavs, nin)}).first;
      }
      const auto& jet_inds = alg->second->Jets();
      Vec4D sum = {0,0,0,0};
      for(const auto& p: moms) sum += p;
      
      double mass2 = 0;
      for(const auto& inds: jet_inds) {
        Vec4D cur = {0,0,0,0};
        for(int ind: inds) cur += moms.at(ind);
        mass2 = std::max(mass2, cur.Abs2());
      }
      return mass2/sum.Abs2();
    }
  };
}

using namespace RESUM;

DECLARE_GETTER(HeavyJetMass,"HeavyJetMass",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,HeavyJetMass>::
operator()(const Parameter_Type &args) const 
{ return new HeavyJetMass(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,HeavyJetMass>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"Heavy Jet Mass"; }

