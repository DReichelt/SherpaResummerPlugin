#include "CParameter.H"


using namespace RESUM;

typedef CParameter_Template<double> CParameter;

DECLARE_GETTER(CParameter,"CParameter",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,CParameter>::
operator()(const Parameter_Type &args) const 
{ return new CParameter(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,CParameter>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"CParameter"; }

