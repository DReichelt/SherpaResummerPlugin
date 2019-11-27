#include "Observables/ColorSinglet/Thrust.H"


using namespace RESUM;

typedef Thrust_Template<double> Thrust2;

DECLARE_GETTER(Thrust2,"Thrust2",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,Thrust2>::
operator()(const Parameter_Type &args) const 
{ return new Thrust2(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,Thrust2>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"Thrust2"; }
