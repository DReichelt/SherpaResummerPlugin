#ifndef RESUM_READER_H
#define RESUM_READER_H
#include <string>
#include <set>
#include <vector>
#include <fstream>
#include <map>
#include "StringTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

namespace RESUM {
  /// Reader class.
  /**
     Class implementing reader to be used for input read in. 
     The first word in every line is interpreted as identifier, everything after it 
     on the same line as value. Multiple words after the first word (as defined by 
     wordseps) are interpreted as vector of values.
  */
  class Reader {
  public:
    /// Separators between words in a line (interpreted as regex).
    std::set<std::string> wordseps = {"[ \t]+"};
    
    /// Separators between individual lines (interpreted as regex, linebreaks are always included).
    std::set<std::string> linseps = {";"};
    /// Comment symbols. Everything after those will be ignored up to the next linebreak.
    std::set<std::string> comments = {"#"};
    
    /// Ignored symbols. These symbols will be stripped from the input.
    std::set<std::string> ignores = {"="};
    
    /// Map identifier -> values.
    std::map<std::string, std::vector<std::string>> values;
    
    /// Path of the input file.
    std::string file_path;
    
    
  public:
    /// \name Standard constructors
    // \{
    /// Default constructor.
    Reader(const std::set<std::string>& wseps={"[ \\t]+"},
           const std::set<std::string>& lseps={";"},
           const std::set<std::string>& comms={"#"},
           const std::set<std::string>& igs={"="}) :
      wordseps(wseps), linseps(lseps), comments(comms), ignores(igs) {}
    
    
    /// Constructor specifying input file. Automatically reads from the file provided.
    Reader(const std::string& file,
           const std::set<std::string>& wordseps={"[ \t]+"},
           const std::set<std::string>& linseps={";"},
           const std::set<std::string>& comments={"#"},
           const std::set<std::string>& ignores={}) :
      file_path(file), wordseps(wordseps), linseps(linseps), comments(comments), ignores(ignores) {
      read();
    }
    // \}
    
    /// Read the file line by line. Has to be called if \p file_path is changed.
    void read();
    
    /// \name Getters.
    // \{
    /// Get a single value, converted to type T.
    /**
       \param name Identifier of the value.
       \param def  Default value, returned if name is not found.
       
       If the identifier \p name points to a vector, the first element will 
       be returned.
    */
    template <typename T>
    T GetValue(const std::string& name, const T& def);
    
    /// Get a vector of values, all converted to type T.
    /**
       \param name Identifier of the vector.
       
       If the identifier is not found, an empty vector is returned.
    */
    template <typename T>
    std::vector<T> GetVector(const std::string& name);
    
    /// Get a matrix of all values, all converted to type T.
    /**
       \param name Identifier of the data.
       \param cols Number of columns in the matrix.
       
       If the identifier is not found, an empty matrix (aka vector<vector<T>>) 
       is returned. The number of rows will be determined by the number of 
       columns and size of the vector in memory. If this size is not a multiple
       cols, the last row will be shorter.
    */
    template <typename T>
    std::vector<std::vector<T>> GetMatrix(const std::string& name, size_t cols);
    // \}
  };

  template <typename T>
  T Reader::GetValue(const std::string& name, const T& def) {
    auto it = values.find(name);
    if(it == values.end() || it->second.empty()) return def;
    return to_type<T>(it->second.at(0));
  }
  

  template <typename T>
  std::vector<T> Reader::GetVector(const std::string& name) {
      auto it = values.find(name);
      if(it == values.end() || it->second.empty()) return std::vector<T>();
      std::vector<T> ret(it->second.size());
      for(size_t i=0; i<ret.size(); i++) ret[i] = to_type<T>(it->second[i]);
      return ret;
  }

  template <typename T>
  std::vector<std::vector<T>> Reader::GetMatrix(const std::string& name, size_t cols) {
    if(cols==0) return std::vector<std::vector<T>>();
    std::vector<T> vals = GetVector<T>(name);
    if(vals.empty()) return std::vector<std::vector<T>>();
    if(vals.size() < cols) return std::vector<std::vector<T>>(1,vals);
    size_t rows = vals.size()/cols;
    std::vector<std::vector<T>> ret(rows);
    for(size_t row=0; row<ret.size()-1; row++) {
      ret[row] = {vals.begin()+cols*row, vals.begin()+cols*(row+1)};
    }
    ret.back() = {vals.begin()+cols*(rows-1), vals.end()};
    return ret;
  }

}
#endif
