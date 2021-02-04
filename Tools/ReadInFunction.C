#include "Tools/ReadInFunction.H"
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



ReadInFunction::ReadInFunction(const std::string& filename, int inc, double yvar) {
  _init(MODE::HERMITE,filename,inc,yvar);
}


ReadInFunction::ReadInFunction(const std::string& filename, double expNLL_LO, 
                     double expNLL_NLO, int inc, double yvar) { 
  m_expNLL_LO = expNLL_LO;
  m_expNLL_NLO = expNLL_NLO;
  _init(MODE::HERMITE,filename,inc,yvar);
}

ReadInFunction::ReadInFunction(const std::string& filename, const std::string& expNLL_LO, 
                               const std::string& expNLL_NLO, int inc, double yvar) {
  m_ID_expNLL_LO = expNLL_LO;
  m_ID_expNLL_NLO = expNLL_NLO;
  _init(MODE::HERMITE,filename,inc,yvar);
}

ReadInFunction::ReadInFunction(const std::string& filename, const std::string& expNLL_LO, 
                               const std::string& expNLL_NLO, const std::string& argFac,
                               int inc, double yvar) {
  m_ID_expNLL_LO = expNLL_LO;
  m_ID_expNLL_NLO = expNLL_NLO;
  m_ID_argFac = argFac;
  _init(MODE::HERMITE,filename,inc,yvar);
}

void ReadInFunction::_init(MODE mode, const std::string& filename, int inc, double yvar) {
  m_mode = mode;
  _read(filename,inc,yvar);
  if(m_mode != MODE::DEFAULT) _calc();
  
}

void ReadInFunction::_read(const std::string& filename, int inc, double yvar) {
  DEBUG_FUNC(filename);
  std::ifstream input(filename);
  if(!input.good()) THROW(fatal_error,"No file " + filename + ".");
  double x;
  double y;
  double yerr;
  std::string row = "";
  getline(input, row);
  if(inc < 1) inc = 1;
  int count = 0;
  while(getline(input, row)){
    const size_t start = row.find_first_not_of(" \t");
    const size_t end = std::min(row.find_first_of("#"), row.find_last_not_of(" \t"));
    if(end <= start) continue;
    if(row.find_first_of("=") == std::string::npos) {
      const std::vector<std::string>& splitrow = split(row.substr(start,end+1),"[ \t]+");
      if(splitrow.size() < 2) THROW(fatal_error,"The file " + filename + " has the wrong format.");
      x = stod(splitrow[0]);
      y = stod(splitrow[1]);
      yerr = splitrow.size() > 2 ? stod(splitrow[2]) : 0;
      
      if(count % inc == 0) {
        m_xvals.push_back(x);
        m_yvals.push_back(y+yvar*yerr);
        m_yerrs.push_back(yerr);
      }
      count++;
    }
    else {
      const std::vector<std::string>& splitrow = split(row.substr(start,end+1),"[ = \t]+");
      if(splitrow[0]==m_ID_expNLL_LO) {
        m_expNLL_LO = to_type<double>(split(splitrow[1]," \t")[0]);
      }
      else if(splitrow[0]==m_ID_expNLL_NLO) {
        m_expNLL_NLO = to_type<double>(split(splitrow[1]," \t")[0]);
      }
      else if(splitrow[0]==m_ID_argFac) {
        m_argFac = to_type<double>(split(splitrow[1]," \t")[0]);
      }
    }
  }
  if(m_xvals.back() != x) {
      m_xvals.push_back(x);
      m_yvals.push_back(y+yvar*yerr);
      m_yerrs.push_back(yerr);
  }
}

double ReadInFunction::operator()(const double x, double& expNLL_NLO) {
  return (*this)(x,m_expNLL_LO,expNLL_NLO);
}

double ReadInFunction::operator()(const double x, double& expNLL_LO, 
                             double& expNLL_NLO) {
  DEBUG_FUNC(x);
  expNLL_LO = m_expNLL_LO*m_argFac;
  expNLL_NLO = m_expNLL_NLO*m_argFac*m_argFac;
  if (m_mode == MODE::DEFAULT) {
    size_t i = 0;
    for(; i<m_xvals.size(); i++) {
      if(m_xvals[i] > x*m_argFac) {
        break;
      }
    }
    if(i == m_xvals.size()) THROW(fatal_error,"No data for requested function. Requested x = "+std::to_string(x)+", x max = "+std::to_string(m_xvals.back()));
    const double x_l = m_xvals[i-1];
    const double y_l = m_yvals[i-1];
    const double y_u = m_yvals[i];
    const double x_u = m_xvals[i];
  
    return y_l+(y_u-y_l)/(x_u-x_l)*(x-x_l);
  }
  else {
    return Interpolate(x*m_argFac); 
  }
}

// void ReadInFunction::PrintYODA(const std::string& filename) {
// #ifdef USING_YODA
//   YODA::Scatter2D plot("/ReadInFunction/ReadInFunction","/ReadInFunction/ReadInFunction");
//   YODA::Scatter2D plot_interp("/ReadInFunction/ReadInFunction","/ReadInFunction/ReadInFunction");
//   double Fvar = stod(ATOOLS::rpa->gen.Variable("RESUM::FFUNCTION::VARIATION","0"));
//   double dummy = 0;  
//   for(size_t i=0; i<m_xvals.size(); i++) {
//     if(ATOOLS::IsZero(m_xvals[i])) continue;
//     plot.addPoint(m_xvals[i],m_yvals[i]-Fvar*m_yerrs[i],0,m_yerrs[i]);
//   }
//   for(size_t i=0; i<m_Rps.size(); i++) {
//     if(ATOOLS::IsZero(m_Rps[i])) continue;
//     plot_interp.addPoint(m_Rps[i],this->operator()(m_Rps[i],dummy));
//   }
  

//   auto maxit = std::max_element(m_xvals.begin(),m_xvals.end());
//   auto minit = std::min_element(m_xvals.begin(),m_xvals.end());
//   for(double x=*minit;x<*maxit;x+=0.01) {
//     if(ATOOLS::IsZero(x)) continue;
//     plot.addPoint(x,0);
//     plot_interp.addPoint(x,this->operator()(x,dummy));
//   }
//   YODA::WriterYODA::write(filename+".yoda",plot);
//   YODA::WriterYODA::write(filename+"_Interpolated_F"+std::to_string(int(Fvar))+".yoda",
//                           plot_interp);
// #endif
// }
