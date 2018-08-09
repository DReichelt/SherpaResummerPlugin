#include "StringTools.H"
#include <fstream>
#include <algorithm>
#include <regex>

using std::set;
using std::string;
using std::vector;

vector<string> RESUM::split(const string& input, const string& regex) {
  // passing -1 as the submatch index parameter performs splitting
  std::regex re(regex);
  std::sregex_token_iterator first{input.begin(), input.end(), re, -1}, last;
  return {first, last};
}

vector<string> RESUM::split(const string& input, const set<string>& vals) {
  string regex = "(";
  for(string v: vals) {
    regex += v;
    regex += "|";
  }
  regex.pop_back();
  regex += ")";
  return split(input, regex);
}
