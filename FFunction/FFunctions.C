#include "FFunction/FFunctions.H"
#include "Tools/Files.H"
#include "ATOOLS/Org/Exception.H"

#include  <fstream>

using namespace RESUM;
using namespace FFUNCTION;



FFunction::FFunction(const std::string& filename) {
  std::ifstream input(FILENAMES::SHARE_DIR+"/FFunction/" + filename);
  if(!input.good()) THROW(fatal_error,"No file " + filename + ".");
  double Rp;
  double F;
  std::string row = "";
  getline(input, row);
  while(getline(input, row)){
    int space1 = row.find(" ");
      int space2 = row.find(" ",space1+1);
      if(space1 == std::string::npos or space2 == std::string::npos) THROW(fatal_error,"The file " + filename + " has wrong format.");
      
      Rp = stod(row.substr(0,space1));
      F = stod(row.substr(space1+1,space2));
      //F_err = row.substr(space2+1);
      Rps.push_back(Rp);
      Fs.push_back(F);          
  }
  input.close(); 
}

double FFunction::operator()(const double Rp) {
  size_t i = 0;
  for(; i<Rps.size(); i++) {
    if(Rps[i] > Rp) {
      break;
    }
  }
  if(i == Rps.size()) THROW(fatal_error,"No data for requested F function. Requested Rp = "+std::to_string(Rp)+", Rp max = "+std::to_string(Rps.back()));
  const double Rp_l = Rps[i-1];
  const double F_l = Fs[i-1];
  const double F_u = Fs[i];
  const double Rp_u = Rps[i];
  return F_l+(F_u-F_l)/(Rp_u-Rp_l)*(Rp-Rp_l);
}

