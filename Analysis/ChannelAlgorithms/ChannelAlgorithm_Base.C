#include "Analysis/ChannelAlgorithms/ChannelAlgorithm_Base.H"

using namespace RESUM;

#define COMPILE__Getter_Function
#define OBJECT_TYPE RESUM::ChAlg
#define PARAMETER_TYPE RESUM::ChAlg_Key
#include "ATOOLS/Org/Getter_Function.C"

#include "ATOOLS/Org/Run_Parameter.H"

std::string ChannelAlgorithm_Base::Channel(const ATOOLS::Particle_List& particles) {
  Vec4D_Vector mom(2+particles.size());
  Flavour_Vector fl(2+particles.size());
  for (size_t i=0; i<particles.size(); i++) {
    mom[2+i]=particles[i]->Momentum();
    fl[2+i]=particles[i]->Flav();
  }
  fl[0]=rpa->gen.Beam1();
  fl[1]=rpa->gen.Beam2();
  return Channel(mom,fl,2);
}
