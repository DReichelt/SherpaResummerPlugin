#include "Observables/TotBroad.H"

using namespace RESUM;
typedef TotBroad_Template<double> TotBroad;

DECLARE_GETTER(TotBroad,"TotBroad",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,TotBroad>::
operator()(const Parameter_Type &args) const 
{ return new TotBroad(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,TotBroad>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"TotBroad"; }

