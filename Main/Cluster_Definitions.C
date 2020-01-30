#include "Main/Cluster_Definitions.H"
#include "ATOOLS/Org/Exception.H"

using namespace RESUM;
using namespace PDF;
using namespace ATOOLS;

Cluster_Definitions::Cluster_Definitions()
{
}

Cluster_Definitions::~Cluster_Definitions()
{
}

CParam Cluster_Definitions::KPerp2
(const Cluster_Amplitude &ampl,int i,int j,int k,
 const Flavour &mo,Mass_Selector *const ms,
 const int kin,const int mode)
{
  THROW(not_implemented,"This function is a dummy.");
  return CParam();
}

Vec4D_Vector Cluster_Definitions::Combine
(const Cluster_Amplitude &ampl,int i,int j,int k,
 const Flavour &mo,Mass_Selector *const ms,
 const int kin,const int mode)
{
  THROW(not_implemented,"This function is a dummy.");
  return Vec4D_Vector();
}
