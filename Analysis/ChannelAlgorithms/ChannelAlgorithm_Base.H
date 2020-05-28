#ifndef RESUM_CHANNELALGORITHM_H
#define RESUM_CHANNELALGORITHM_H
#include <vector>
#include <string>
#include <memory>
#include "Tools/StringTools.H"
#include "Tools/Key_Base.H"

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Org/Message.H"


namespace RESUM {

  // struct ChannelAlgorithm_Key {
  //   ChannelAlgorithm_Key(const std::string& key)
  //     : m_key(key) {}
  //   ChannelAlgorithm_Key(const std::string& key,
  //                        std::vector<std::string> params)
  //     : m_key(key), m_params(params) {
  //     // conventionally remove key from parmas if present
  //     if(m_params[0]==m_key) m_params.erase(m_params.begin());
  //   }
  //   std::string m_key = "ChannelAlgorithm";
  //   std::vector<std::string> m_params;
  // };
  
  // convenience typedef
  typedef Key_Base ChAlg_Key;

  class ChannelAlgorithm_Base {
  public:
    typedef std::shared_ptr<ChannelAlgorithm_Base> Ptr;
    ChannelAlgorithm_Base() = default;
    ChannelAlgorithm_Base(const ChAlg_Key& params)
      : m_name(params.Name()),  m_params(params.m_params) { }
  
    virtual std::string Channel(const std::vector<ATOOLS::Vec4D>& ip,
                                const std::vector<ATOOLS::Flavour>& fl,
                                const size_t &nin) = 0;

    std::string Channel(const ATOOLS::Particle_List& particles);
    
    virtual std::vector<std::string> ChannelNames() = 0;
    
  protected:
    std::string m_name;
    std::vector<std::string> m_params;

  }; // end of class ChannelAlgorithm_Base

  // convenience typedef
  typedef ChannelAlgorithm_Base ChAlg;

  // prepare for Getter implementation
  typedef ATOOLS::Getter_Function <ChAlg,ChAlg_Key> ChAlg_Getter;
  

} // end of namespace RESUM

#endif
