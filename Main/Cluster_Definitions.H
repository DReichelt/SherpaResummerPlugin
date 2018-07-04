#ifndef RESUM__Main__Cluster_Definitons_H
#define RESUM__Main__Cluster_Definitons_H

#include "PDF/Main/Cluster_Definitions_Base.H"

namespace RESUM {

  class Cluster_Definitions: public PDF::Cluster_Definitions_Base {
  public:

    Cluster_Definitions();

    ~Cluster_Definitions();

    PDF::CParam KPerp2
    (const ATOOLS::Cluster_Amplitude &ampl,int i,int j=-1,int k=-1,
     const ATOOLS::Flavour &mo=ATOOLS::Flavour(kf_none),
     ATOOLS::Mass_Selector *const ms=NULL,
     const int kin=-1,const int mode=0);

    ATOOLS::Vec4D_Vector Combine
    (const ATOOLS::Cluster_Amplitude &ampl,int i,int j,int k=-1,
     const ATOOLS::Flavour &mo=ATOOLS::Flavour(kf_none),
     ATOOLS::Mass_Selector *const ms=NULL,
     const int kin=0,const int mode=0);

    PDF::CParam ResumScale(ATOOLS::Cluster_Amplitude *const ampl);

  };// end of class Cluster_Definitions

}// end of namespace RESUM

#endif
