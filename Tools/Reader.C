#include "Reader.H"

using std::vector;
using std::string;

namespace RESUM {
  void Reader::read() {
    msg_Debugging()<<"Read file "<<file_path<<"\n";
    string cmms;
    for(const string& cm: comments) cmms += cm;
    std::ifstream input(file_path);
    string row;
    if(!input.good()) THROW(fatal_error,"No file "+file_path+".");
    while(getline(input, row)) {
      row = row.substr(0,row.find_first_of(cmms));
      for(const string& ign: ignores) {
        for(size_t pos=row.find(ign); pos!=std::string::npos; pos=row.find(ign,pos)) {
          row.erase(pos,ign.size());
        }
      }
      if(!row.empty()) {
        vector<string> lines = split(row,linseps);
        for(const string& line: lines) {
          vector<string> words = split(line,wordseps);
          if(!words.empty()) {
            values[words[0]] = {words.begin()+1, words.end()};
          }
        }
      }
    }
  }
}
