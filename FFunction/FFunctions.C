#include "FFunction/FFunctions.H"
#include "Tools/Files.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/MathTools.H"



#include  <fstream>

using namespace RESUM;
using namespace FFUNCTION;



FFunction::FFunction(const std::string& filename) {
  DEBUG_FUNC(filename);
  std::ifstream input(FILENAMES::SHARE_DIR+"/FFunction/" + filename);
  if(!input.good()) THROW(fatal_error,"No file " + filename + ".");
  m_mode = MODE::HERMITE;
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

      // double w = 1.;
      // for(size_t i=0; i<m_bweights.size(); i++) {
      //   const double wi = m_Rps[i]-Rp;
      //   if(ATOOLS::IsZero(wi)) THROW(fatal_error, "Points in FFunction must not be equal!")
      //   m_bweights[i] /= wi;
      //   w /= -wi;
      // }
      // m_bweights.push_back(w);
      m_xvals.push_back(Rp);
      m_yvals.push_back(F);          
  }
  if(m_mode != MODE::DEFAULT) _calc();
  input.close(); 
}

double FFunction::operator()(const double Rp, double& FexpNLL_NLO) {
  DEBUG_FUNC(Rp);
  if (m_mode == MODE::DEFAULT) {
    size_t i = 0;
    for(; i<m_xvals.size(); i++) {
      if(m_xvals[i] > Rp) {
        break;
      }
    }
    if(i == m_xvals.size()) THROW(fatal_error,"No data for requested F function. Requested Rp = "+std::to_string(Rp)+", Rp max = "+std::to_string(m_xvals.back()));
    const double Rp_l = m_xvals[i-1];
    const double F_l = m_yvals[i-1];
    const double F_u = m_yvals[i];
    const double Rp_u = m_xvals[i];
  
    return F_l+(F_u-F_l)/(Rp_u-Rp_l)*(Rp-Rp_l);
  }
  else {
    return Interpolate(Rp); 
  }
}

// double FFunction::LaplacePol(const double Rp) const {
//   double num = 0;
//   double den = 0;
//   for(size_t i=0; i<m_xvals.size(); i++) {
//     if(ATOOLS::IsZero(Rp-m_xvals[i],1e-6)) {
//       return m_yvals[i];
//     }
//     const double w = m_bweights[i]/(Rp-m_xvals[i]);
//     num += w*m_yvals[i];
//     den += w;
//   }
//   return num/den;
// }
