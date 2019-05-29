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
        m_CF = m_NC;
        m_CA = m_NC;
        msg_Debugging()<<"Calculate QCD parameters in large NC limit.\n";
      }
      else msg_Debugging()<<"Calculate QCD parameters with NC = 3.\n";
    }
    double NC() const {return m_NC;}
    double CA() const {return m_CA;}
    double CF() const {return m_CF;}
    double TR() const {return m_TR;}
    double beta0(double scale2) const;
    double beta1(double scale2) const;
    double K_CMW(double scale2) const;
    double CollDimGlue(double scale2) const;
    double CollDimQuark(double scale2) const;
  };
  
  const double s_Nc = 3.;
  const double s_CF = (s_Nc*s_Nc-1.)/(2.*s_Nc);
  const double s_CA = s_Nc;
  const double s_TR = 1./2.;
  const double s_eps = .000001;
  
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
    std::function<double(double)> m_F;
    std::vector<int> signlabels;
    int n_g, n_q, n_aq, color_sings;
    double s_12, m_rn[2];

    /// Leg ids ordered as quarks-antiquarks-gluons in n-parton amplitude.
    std::vector<size_t> m_ordered_ids;
    
    enum MODE {
               RESUM = 0,
               PDFEXPAND = 1  << 0,
               COLLEXPAND = 1 << 1,
               SOFTEXPAND = 1 << 2,
               LARGENC = 1    << 3,
               EXPAND = COLLEXPAND | SOFTEXPAND
    };

    const std::map<std::string,MODE> m_ModeToEnum = {{"RESUM", MODE::RESUM},
                                                     {"PDFEXPAND", MODE::PDFEXPAND},
                                                     {"COLLEXPAND", MODE::COLLEXPAND},
                                                     {"SOFTEXPAND", MODE::SOFTEXPAND},
                                                     {"LARGENC", MODE::LARGENC},
                                                     {"EXPAND", MODE::EXPAND}};
    
    MODE m_amode = MODE::RESUM;
    
    double CalcS(const double L, const double LResum, double &Softexp);
    double CalcColl(const double L, const double LResum, const int order,
		    double &Rp,double &CollexpLL, double& CollexpNLL);
    double CalcPDF(const double L, const double LResum, double &PDFexp);
    double T(const double x);

    double Value(const double v, const double LResum, double epRatio);

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

    std::pair<int, double> Result(int n) const { return m_ress.at(n); }
    virtual void ResetObservables();

  };// end of class Resum
  
}// end of namespace RESUM

#endif
