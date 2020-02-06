#include "FFunction/FFunctions.H"
#include "Tools/Files.H"
#include "Tools/StringTools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Run_Parameter.H"

#ifdef USING_YODA
#include "YODA/Scatter2D.h"
#include "YODA/WriterYODA.h"
#endif

#include  <fstream>

using namespace RESUM;
using namespace FFUNCTION;



FFunction::FFunction(const std::string& filename) {
  m_mode = MODE::HERMITE;
  std::string name = filename;
  std::ifstream input(FILENAMES::SHARE_DIR+"/FFunction/" + name);
  if(!input.good()) THROW(fatal_error,"No file " + name + ".");
  std::string line = "";
  while(getline(input,line)) {
    if(line.find("GOTO") == 0) {
      name = line.substr(line.find(" ")+1);
      input.close();
      input = std::ifstream(FILENAMES::SHARE_DIR+"/FFunction/" + name);
    }
  }
  int inc = std::stoi(ATOOLS::rpa->gen.Variable("RESUM::FFUNCTION::INC","1"));
  _read(name,inc);
  if(m_mode != MODE::DEFAULT) _calc();
  input.close(); 
  if(ATOOLS::rpa->gen.Variable("RESUM::FFUNCTION::PLOT","0") != "0") {
    size_t i = filename.rfind(".",filename.size());
    PrintYODA(filename.substr(0,i)+"_Inc"+std::to_string(inc));
  }
}


FFunction::FFunction(const std::string& filename, double F2) : FFunction(filename) {
  m_F2 = F2;
}

void FFunction::_read(const std::string& filename,int inc) {
  DEBUG_FUNC(filename);
  std::ifstream input(FILENAMES::SHARE_DIR+"/FFunction/" + filename);
  if(!input.good()) THROW(fatal_error,"No file " + filename + ".");
  double Fvar = stod(ATOOLS::rpa->gen.Variable("RESUM::FFUNCTION::VARIATION","0"));

  double Rp;
  double F;
  double Ferr;
  std::string row = "";
  getline(input, row);
  if(inc < 1) inc = 1;
  int count = 0;
  while(getline(input, row)){
    const std::vector<std::string>& splitrow = split(row,"[ \t]+");
    if(splitrow.size() < 2) THROW(fatal_error,"The file " + filename + " has the wrong format.");
    Rp = stod(splitrow[0]);
    F = stod(splitrow[1]);
    Ferr = splitrow.size() > 2 ? stod(splitrow[2]) : 0;
      
    m_Rps.push_back(Rp);

    if(count % inc == 0) {

      // double w = 1.;
      // for(size_t i=0; i<m_bweights.size(); i++) {
      //   const double wi = m_Rps[i]-Rp;
      //   if(ATOOLS::IsZero(wi)) THROW(fatal_error, "Points in FFunction must not be equal!")
      //   m_bweights[i] /= wi;
      //   w /= -wi;
      // }
      // m_bweights.push_back(w);
      m_xvals.push_back(Rp);
      m_yvals.push_back(F+Fvar*Ferr);
      m_yerrs.push_back(Ferr);
    }
    count++;
  }
  if(m_xvals.back() != Rp) {
      m_xvals.push_back(Rp);
      m_yvals.push_back(F+Fvar*Ferr);
      m_yerrs.push_back(Ferr);
  }
}

double FFunction::operator()(const double Rp, double& FexpNLL_NLO) {
  DEBUG_FUNC(Rp);
  FexpNLL_NLO = m_F2;
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

void FFunction::PrintYODA(const std::string& filename) {
#ifdef USING_YODA
  YODA::Scatter2D plot("/FFunction/FFunction","/FFunction/FFunction");
  YODA::Scatter2D plot_interp("/FFunction/FFunction","/FFunction/FFunction");
  double Fvar = stod(ATOOLS::rpa->gen.Variable("RESUM::FFUNCTION::VARIATION","0"));
  double dummy = 0;  
  for(size_t i=0; i<m_xvals.size(); i++) {
    if(ATOOLS::IsZero(m_xvals[i])) continue;
    plot.addPoint(m_xvals[i],m_yvals[i]-Fvar*m_yerrs[i],0,m_yerrs[i]);
  }
  for(size_t i=0; i<m_Rps.size(); i++) {
    if(ATOOLS::IsZero(m_Rps[i])) continue;
    plot_interp.addPoint(m_Rps[i],this->operator()(m_Rps[i],dummy));
  }
  

  auto maxit = std::max_element(m_xvals.begin(),m_xvals.end());
  auto minit = std::min_element(m_xvals.begin(),m_xvals.end());
  for(double x=*minit;x<*maxit;x+=0.01) {
    if(ATOOLS::IsZero(x)) continue;
    plot.addPoint(x,0);
    plot_interp.addPoint(x,this->operator()(x,dummy));
  }
  YODA::WriterYODA::write(filename+".yoda",plot);
  YODA::WriterYODA::write(filename+"_Interpolated_F"+std::to_string(int(Fvar))+".yoda",
                          plot_interp);
#endif
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
