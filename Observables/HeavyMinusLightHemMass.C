#include "Observables/HeavyMinusLightHemMass.H"
using namespace RESUM;
typedef HeavyMinLightMass_Template<double> HeavyMinLightMass;

DECLARE_GETTER(HeavyMinLightMass,"HeavyMinLightMass",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,HeavyMinLightMass>::
operator()(const Parameter_Type &args) const 
{ return new HeavyMinLightMass(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,HeavyMinLightMass>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"HeavyMinLightMass"; }

