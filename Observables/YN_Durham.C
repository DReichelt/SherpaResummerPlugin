#include "Observables/YN_Durham.H"
using namespace RESUM;

typedef YN_Durham<3> Y3_Durham;
DECLARE_GETTER(Y3_Durham,"Y3_Durham",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,Y3_Durham>::
operator()(const Parameter_Type &args) const 
{ return new Y3_Durham(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,Y3_Durham>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"Y3_Durham"; }


typedef YN_Durham<4> Y4_Durham;
DECLARE_GETTER(Y4_Durham,"Y4_Durham",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,Y4_Durham>::
operator()(const Parameter_Type &args) const 
{ return new Y4_Durham(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,Y4_Durham>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"Y4_Durham"; }


typedef YN_Durham<5> Y5_Durham;
DECLARE_GETTER(Y5_Durham,"Y5_Durham",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,Y5_Durham>::
operator()(const Parameter_Type &args) const 
{ return new Y5_Durham(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,Y5_Durham>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"Y5_Durham"; }

typedef YN_Durham<6> Y6_Durham;
DECLARE_GETTER(Y6_Durham,"Y6_Durham",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,Y6_Durham>::
operator()(const Parameter_Type &args) const 
{ return new Y6_Durham(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,Y6_Durham>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"Y6_Durham"; }

typedef YN_Durham<7> Y7_Durham;
DECLARE_GETTER(Y7_Durham,"Y7_Durham",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,Y7_Durham>::
operator()(const Parameter_Type &args) const 
{ return new Y7_Durham(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,Y7_Durham>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"Y7_Durham"; }

typedef YN_Durham<8> Y8_Durham;
DECLARE_GETTER(Y8_Durham,"Y8_Durham",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,Y8_Durham>::
operator()(const Parameter_Type &args) const 
{ return new Y8_Durham(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,Y8_Durham>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"Y8_Durham"; }
