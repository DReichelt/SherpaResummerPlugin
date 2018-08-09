#ifndef RESUM_STRING_TOOLS
#define RESUM_STRING_TOOLS
/** \file 
    Some basic functions involving strings.
 */
#include <string>
#include <set>
#include <vector>
#include <sstream>

namespace RESUM {
  /// Split string \p input into multiple strings
  /**
     \param input The string to be split.
     \param regex Regular expression, to determine where to split string.
     
     The string will be split at every point that matches the regex.
  */
  std::vector<std::string> split(const std::string& input, const std::string& regex);

  /// Split string \p input into multiple strings
  /**
     \param input The string to be split.
     \param vals  Container with values of seperators.
     
     The string will be split every time a seperator is found.
  */
  std::vector<std::string> split(const std::string& input, const std::set<std::string>& vals);
  
  /// Convert string to type T.
  /**
     \param value String to be converted.
  */
  template <typename T>
  T to_type(const std::string& value)  {
    if(std::is_same<T,int>::value) return std::stoi(value);
    if(std::is_same<T,double>::value) return std::stod(value);
    if(std::is_same<T,long>::value) return std::stol(value);
    if(std::is_same<T,long long>::value) return std::stoll(value);
    if(std::is_same<T,float>::value) return std::stof(value);
    if(std::is_same<T,long double>::value) return std::stold(value);
    
    std::istringstream conv(value);
    T ret;
    conv>>ret;
    return ret;
  }   
}
#endif
