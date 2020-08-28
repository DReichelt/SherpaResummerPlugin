#include "Analysis/ChannelAlgorithms/ChannelAlgorithm_Base.H"

using namespace RESUM;

#define COMPILE__Getter_Function
#define OBJECT_TYPE RESUM::ChAlg
#define PARAMETER_TYPE RESUM::ChAlg_Key
#include "ATOOLS/Org/Getter_Function.C"

#include "ATOOLS/Org/Run_Parameter.H"

std::string ChannelAlgorithm_Base::Channel(const ATOOLS::Particle_List& particles, bool addTag, 
                                           std::vector<ATOOLS::Vec4D>* pout) {
  Vec4D_Vector mom(2+particles.size());
  Flavour_Vector fl(2+particles.size());
  for (size_t i=0; i<particles.size(); i++) {
    mom[2+i]=particles[i]->Momentum();
    fl[2+i]=particles[i]->Flav();
  }
  fl[0]=rpa->gen.Beam1();
  fl[1]=rpa->gen.Beam2();
  return Channel(mom,fl,2,addTag,pout);
}

std::string ChannelAlgorithm_Base::Channel(const std::vector<ATOOLS::Vec4D>& ip,
                                           const std::vector<ATOOLS::Flavour>& fl,
                                           const size_t &nin, bool addTag,
                                           std::vector<ATOOLS::Vec4D>* pout) {
  std::string ch = Channel(ip,fl,nin,pout);
  if(addTag and m_tag != "") ch += "_"+m_tag;
  return ch;
}

std::vector<std::string> ChannelAlgorithm_Base::ChannelNames(bool addTag) {
  std::vector<std::string> chs = ChannelNames();
  if(addTag and m_tag != "") {
    for(std::string& ch: chs) ch += "_"+m_tag;
  }
  return chs;
}
