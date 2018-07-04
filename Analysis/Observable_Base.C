#include "Analysis/Observable_Base.H"

using namespace RESUM;

#define COMPILE__Getter_Function
#define OBJECT_TYPE RESUM::Observable_Base
#define PARAMETER_TYPE RESUM::Observable_Key
#include "ATOOLS/Org/Getter_Function.C"

#include "ATOOLS/Math/Poincare.H"

Observable_Base::Observable_Base(const Observable_Key &args):
  m_name(args.m_name)
{
}

Observable_Base::~Observable_Base()
{
}

double Observable_Base::Value(ATOOLS::Particle_List *const pl)
{
  std::vector<Vec4D> p(pl->size());
  std::vector<Flavour> fl(pl->size());
  for (size_t i(0);i<pl->size();++i) {
    p[i]=(*pl)[i]->Momentum();
    fl[i]=(*pl)[i]->Flav();
  }
  return Value(&p.front(),&fl.front(),pl->size(),0);
}

double Observable_Base::Shift(ATOOLS::NLO_subevt *sub)
{
  Vec4D p(sub->p_mom[sub->m_ijt]);
  Poincare cms(sub->p_mom[0]+sub->p_mom[1]);
  cms.Boost(p);
  Obs_Params ps(Parameters(sub->p_mom,sub->p_fl,sub->m_n,sub->m_ijt));
  return ps.m_d-ps.m_b/2.0*log(sqr(2.0*p[0])/sub->m_mu2[stp::res]);
}
