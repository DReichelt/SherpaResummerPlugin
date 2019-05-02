#include <algorithm>
using std::max;
#include "ThrustMinor.H"


using namespace RESUM;

typedef ThrustMinor_Template<double> ThrustMinor;

DECLARE_GETTER(ThrustMinor,"ThrustMinor",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,ThrustMinor>::
operator()(const Parameter_Type &args) const 
{ return new ThrustMinor(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,ThrustMinor>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"Thrust Major"; }

