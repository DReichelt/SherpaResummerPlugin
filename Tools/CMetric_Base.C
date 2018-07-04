#include "Tools/CMetric_Base.H"

using namespace RESUM;

#define COMPILE__Getter_Function
#define OBJECT_TYPE RESUM::CMetric_Base
#define PARAMETER_TYPE RESUM::CMetric_Key
#include "ATOOLS/Org/Getter_Function.C"


CMetric_Base::CMetric_Base()
{
}

CMetric_Base::CMetric_Base(const CMetric_Key &args): m_name(args.m_name)
{
}

void CMetric_Base::CalcMetric() 
{
}


void CMetric_Base::CalcIMetric() 
{
}

void CMetric_Base::CalcTs() 
{
}


CMetric_Base::~CMetric_Base()
{
}



typedef ATOOLS::Getter_Function
<CMetric_Base,CMetric_Key> CMetric_Getter;

CMetric_Base* CMetric_Base::GetCM(const CMetric_Key &args)
{
  const CMetric_Getter::Getter_List &glist(CMetric_Getter::GetGetters());
  for (CMetric_Getter::Getter_List::const_iterator git(glist.begin());
       git!=glist.end();++git) {
    CMetric_Base *cm=(*git)->GetObject(args);
    if (cm) return cm;
  }
  return NULL;
}

