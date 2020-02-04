#include "Observables/WideBroad.H"
using namespace RESUM;
typedef WideBroad_Template<double> WideBroad;

DECLARE_GETTER(WideBroad,"WideBroad",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,WideBroad>::
operator()(const Parameter_Type &args) const 
{ return new WideBroad(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,WideBroad>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"WideBroad"; }

