#include "Tools/Key_Base.H"
#include "Tools/StringTools.H"

using namespace RESUM;

Key_Base::Key_Base(const std::string& name) : m_name(name) {}


Key_Base::Key_Base(const std::string& name,
                   std::vector<std::string> params)
  : m_name(name), m_params(params) {
  // conventionally remove key from parmas if present
  if(m_params[0]==m_name) m_params.erase(m_params.begin());
  for(const std::string& p: m_params) {
    std::vector<std::string> keyVal = split(p,":");
    if(keyVal.size() > 1) {
      m_params_map[keyVal[0]] = keyVal[1];
    }
  }
}


