#include "Observables/Oblateness.H"

using namespace RESUM;
typedef Oblateness_Template<double> Oblateness;

DECLARE_GETTER(Oblateness,"Oblateness",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,Oblateness>::
operator()(const Parameter_Type &args) const 
{ return new Oblateness(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,Oblateness>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"Oblateness"; }

