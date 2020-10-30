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
  : m_name(name) {
  DEBUG_FUNC(name);
  // conventionally remove key from parmas if present
  // m_params.clear();
  // m_params.reserve(params.size());
  if(params[0]==m_name) params.erase(params.begin());
  for(std::string p: std::move(params)) {
    msg_Debugging()<<p<<"\n";
    for(std::string pp: split(p,std::set<std::string>({" ",","}))) { 
      msg_Debugging()<<pp<<"\n";
      m_params.emplace_back(std::move(pp));
      std::vector<std::string> keyVal = split(m_params.back(),std::set<std::string>({":"}));
      if(keyVal.size() > 1) {
        msg_Debugging()<<keyVal[0]<<" -> "<<keyVal[1]<<"\n";
        m_params_map[keyVal[0]] = keyVal[1];
      }
      
    }
  }
}



