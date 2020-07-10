#include "Tools/Key_Base.H"
#include "Tools/StringTools.H"
#include "ATOOLS/Org/Message.H"

using namespace RESUM;

Key_Base::Key_Base(const std::string& name) : m_name(name) {
  DEBUG_FUNC(name);
  msg_Debugging()<<"No parameters.";
}


Key_Base::Key_Base(const std::string& name,
                   std::vector<std::string> params)
  : m_name(name), m_params(params) {
  DEBUG_FUNC(name);
  // conventionally remove key from parmas if present
  if(m_params[0]==m_name) m_params.erase(m_params.begin());
  for(const std::string& p: m_params) {
    msg_Debugging()<<p<<"\n";
    std::vector<std::string> keyVal = split(p,std::set<std::string>({":"," ",","}));
    if(keyVal.size() > 1) {
      msg_Debugging()<<keyVal[0]<<" -> "<<keyVal[1]<<"\n";
      m_params_map[keyVal[0]] = keyVal[1];
    }
  }
}



