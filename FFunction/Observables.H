#ifndef RESUM_FFUNCOBSERVABLES_H
#define RESUM_FFUNCOBSERVABLES_H


#include "MultiplePrecision.H"
#include "Analysis/Observable_Base.H"

using mp::max;
#include "Observables/HeavyJetMass.H"




namespace RESUM {
  typedef Observable_Base_Template<mpfloat> Observable_Base_MP;


  typedef HeavyJetMass_Template<mpfloat> HeavyJetMassMP;

  typedef ATOOLS::Getter_Function
    <Observable_Base_MP,Observable_Key> Observable_Getter_MP;

}

using namespace RESUM;

#endif
