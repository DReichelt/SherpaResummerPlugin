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
    void SetKwArg(const std::string& kw, const std::string& arg) {m_params_map[kw] = arg;}
    void SetKwArg(const std::map<std::string,std::string>& kw_args) {for(auto& kw_arg: kw_args) SetKwArg(kw_arg.first, 
                                                                                                         kw_arg.second);}
    void SetKwArg(const Key_Base& kw_args) {SetKwArg(kw_args.m_params_map);}
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
