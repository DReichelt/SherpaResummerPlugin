#include "FFunction/Observables.H"

using namespace RESUM;

#define COMPILE__Getter_Function
#define OBJECT_TYPE Observable_Base_MP
#define PARAMETER_TYPE Observable_Key
#include "ATOOLS/Org/Getter_Function.C"


DECLARE_GETTER(HeavyJetMassMP,"HeavyJetMass",Observable_Base_MP,Observable_Key);
Observable_Base_MP *ATOOLS::Getter<Observable_Base_MP,Observable_Key,HeavyJetMassMP>::
  operator()(const Parameter_Type &args) const 
{ return new HeavyJetMassMP(args); }
void ATOOLS::Getter<Observable_Base_MP,Observable_Key,HeavyJetMassMP>::
  PrintInfo(std::ostream &str,const size_t width) const
{ str<<"Heavy Jet Mass"; }

