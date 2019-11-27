#include "Observables/EParameter.H"

using namespace RESUM;
typedef EParameter_Template<double> EParameter;

DECLARE_GETTER(EParameter,"EParameter",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,EParameter>::
operator()(const Parameter_Type &args) const 
{ return new EParameter(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,EParameter>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"EParameter"; }
