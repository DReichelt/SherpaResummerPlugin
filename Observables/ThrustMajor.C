#include <algorithm>
using std::max;
#include "ThrustMajor.H"


using namespace RESUM;

typedef ThrustMajor_Template<double> ThrustMajor;

DECLARE_GETTER(ThrustMajor,"ThrustMajor",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,ThrustMajor>::
operator()(const Parameter_Type &args) const 
{ return new ThrustMajor(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,ThrustMajor>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"Thrust Major"; }

