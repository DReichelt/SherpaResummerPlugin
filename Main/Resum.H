#ifndef RESUM_RESUM_H
#define RESUM_RESUM_H

#include "PDF/Main/Shower_Base.H"

#include "Main/Comix_Interface.H"
#include "Tools/CMetric_Base.H"

#include <map>
#include <string>
#include <functional>

namespace MODEL {
  class Running_AlphaS;
}

namespace PDF {
  class PDF_Base;
}


namespace RESUM {


  struct Params {
  private:
    double m_NC = 3;
    double m_CF = (m_NC*m_NC-1.)/(2.*m_NC);
    double m_CA = m_NC;
    double m_TR = 1./2.;
    bool m_largeNC = false;
    MODEL::Running_AlphaS *p_as;
  public:
    Params() = default;
    Params(MODEL::Running_AlphaS *as, bool largeNC=false) : m_largeNC(largeNC) {
      p_as = as;
      if(m_largeNC) {
        m_CF = m_TR*m_NC;
        m_CA = m_NC;
        msg_Debugging()<<"Calculate QCD parameters in large NC limit.\n";
      }
      else msg_Debugging()<<"Calculate QCD parameters with NC = 3.\n";
    }
    double NC() const {return m_NC;}
    double CA() const {return m_CA;}
    double CF() const {return m_CF;}
    double TR() const {return m_TR;}
    bool largeNC() const {return m_largeNC;}
    double beta0(double scale2) const;
    double beta1(double scale2) const;
    double K_CMW(double scale2) const;
    double CollDimGlue(double scale2) const;
    double CollDimQuark(double scale2) const;
  };
    
  class Cluster_Definitions;
  template <typename T>
  class Observable_Base_Template;
  typedef Observable_Base_Template<double> Observable_Base;
  
  class Resum: public PDF::Shower_Base {
  private:

    Comix_Interface m_comix;

    Cluster_Definitions *p_clus;

    ATOOLS::Cluster_Amplitude *p_ampl;

    PDF::PDF_Base    ** p_pdf;
    
    MODEL::Running_AlphaS *p_as;

    Params m_params;

    //list of available and initialized bases
    CMetric_Map    m_cmetrics;
    //suitable basis for actual event
    CMetric_Base * p_cmetric;
    
    std::vector<Observable_Base*>     m_obss;
    std::vector<ATOOLS::Histogram*>   m_hist;
    std::vector<std::pair<double,double> > m_ress;

    std::vector<ATOOLS::Vec4D> m_cms, momlabels;
    std::vector<ATOOLS::Flavour> flavlabels;
    std::vector<double> m_Qij, m_a, m_b, m_logdbar;
    std::function<double(double, double&)> m_F;
    std::vector<int> signlabels;
    int n_g, n_q, n_aq, color_sings;
    double s_12, m_rn[2];
    size_t m_n;
    bool m_init = false;
    
    /// Leg ids ordered as quarks-antiquarks-gluons in n-parton amplitude.
    std::vector<size_t> m_ordered_ids;
    
    enum MODE {
               RESUM = 0,               // calculate resummed cumulative
               PDFEXPAND = 1  << 0,     // include pdf in expansion
               COLLEXPAND = 1 << 1,     // include collinear pieces in expansion
               SOFTEXPAND = 1 << 2,     // include soft pieces in expansion
               LARGENC = 1    << 3,     // calculate in large nc limit
               CKINV = 1      << 4,     // check the inverse of metric
               CKCOUL = 1     << 5,     // check effect of coulomb phases
               EXPAND = COLLEXPAND | SOFTEXPAND | PDFEXPAND // conveninece def.
    };

    const std::map<std::string,MODE> m_ModeToEnum = {{"RESUM", MODE::RESUM},
                                                     {"PDFEXPAND", MODE::PDFEXPAND},
                                                     {"COLLEXPAND", MODE::COLLEXPAND},
                                                     {"SOFTEXPAND", MODE::SOFTEXPAND},
                                                     {"LARGENC", MODE::LARGENC},
                                                     {"EXPAND", MODE::EXPAND},
                                                     {"CKINV", MODE::CKINV},
                                                     {"CKCOUL", MODE::CKCOUL},
                                                     {"CHECKINVERSE", MODE::CKINV},
                                                     {"CHECKCOULOMB", MODE::CKCOUL}};

    enum MATCH_MODE {
      NONE = 0,
      LO   = 1          << 1,
      NLO = 1       << 2,
      ADD  = 1          << 3,
    };

    const std::map<std::string,MATCH_MODE> m_MModeToEnum = {{"NONE", MATCH_MODE::NONE},
                                                     {"LO", MATCH_MODE::LO},
                                                     {"NLO", MATCH_MODE::NLO},
                                                     {"ADD", MATCH_MODE::ADD}};

      
    MODE m_amode = MODE::RESUM;
    MATCH_MODE m_mmode = MATCH_MODE::NONE;
    
    double CalcS(const double L, const double LResum, double &SoftexpNLL_LO, double& SoftexpNLL_NLO, MODE Check=MODE::RESUM);
    double CalcColl(const double L, const double LResum, const int order,
		    double &Rp, MatrixD& G, double& S1);
    double CalcPDF(const double L, const double LResum, double &PDFexp);
    double T(const double x);

    double Value(const double v, const double LResum, double epRatio);
    void FillValue(size_t i, const double v, const double LResum, double epRatio);
    
  public:

    Resum(PDF::ISR_Handler *const isr,
	  MODEL::Model_Base *const model);

    ~Resum();


    int  PerformShowers();
    int  PerformDecayShowers();

    bool ExtractPartons(ATOOLS::Blob_List *const bl);
    void CleanUp();

    PDF::Cluster_Definitions_Base *GetClusterDefinitions();

    bool PrepareShower(ATOOLS::Cluster_Amplitude *const ampl,
		       const bool & soft=false);

   
    double CplFac(const ATOOLS::Flavour &fli,const ATOOLS::Flavour &flj,
		  const ATOOLS::Flavour &flk,const int type,
		  const int cpl,const double &mu2) const;

    virtual size_t AddObservable(Observable_Base *const obs,
			 ATOOLS::Histogram *const h);

    virtual size_t AddObservable(const std::string& name,
                                 const std::vector<double>& xvals);
   
    size_t nLegs() const {return n_g+n_q+n_aq+color_sings;}
    std::pair<int, double> Result(int n) const { return m_ress.at(n); }
    bool Initialized() const { return m_init; } 
    virtual void ResetObservables();

    std::vector<std::vector<double>> m_resNLL;
    std::vector<std::vector<double>> m_resExpLO;
    std::vector<std::vector<double>> m_resExpNLO;
    std::vector<std::vector<double>> m_xvals;
    
  };// end of class Resum
  
}// end of namespace RESUM

#endif
