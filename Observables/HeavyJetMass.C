#include <algorithm>
using std::max;
#include "HeavyJetMass.H"


using namespace RESUM;

typedef HeavyJetMass_Template<double> HeavyJetMass;

DECLARE_GETTER(HeavyJetMass,"HeavyJetMass",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,HeavyJetMass>::
operator()(const Parameter_Type &args) const 
{ return new HeavyJetMass(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,HeavyJetMass>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"Heavy Jet Mass"; }

