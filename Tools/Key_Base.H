#ifndef RESUM_KEYBASE_H
#define RESUM_KEYBASE_H

#include <string>
#include <vector>
#include <map>

namespace RESUM {

  struct Key_Base {
    std::vector<std::string> m_params;
    std::map<std::string,std::string> m_params_map;
  public:
    Key_Base(const std::string& name);
    Key_Base(const std::string& name,
             std::vector<std::string> params);
    std::string m_name;
    const std::string& Name() const {return m_name;}
    const std::string& Arg(size_t i) const {return m_params.at(i);}
    const std::string& KwArg(const std::string& kw) const {return m_params_map.at(kw);}
    const std::string& KwArg(const std::string& kw,
                             const std::string& def) const {
      auto ret = m_params_map.find(kw);
      if(ret == m_params_map.end()) {
        return def;
      }
      else {
        return ret->second;
      }
    }
  };

}

#endif
