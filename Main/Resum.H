#ifndef RESUM_RESUM_H
#define RESUM_RESUM_H

#include "PDF/Main/Shower_Base.H"

#include "Main/Comix_Interface.H"
#include "Main/Params.H"
#include "Tools/CMetric_Base.H"
#include "Analysis/Observable_Base.H"

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

    PDF::ISR_Handler* p_isr;

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
    std::vector<double> m_Qij, m_a, m_b, m_logdbar, m_etamin;
    std::vector<std::pair<size_t, size_t>> m_kij;
    std::vector<int> m_deltad;
    std::function<double(double, double&)> m_F;
    std::function<double(double, double&)> m_Sngl;
    std::vector<int> signlabels;
    int n_g, n_q, n_aq, color_sings;
    double s_12, m_rn[2];

    // grooming parameters
    double m_zcut, m_beta;
    
    double m_logpow;

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
               IGNCOLL = 1    << 6,     // ignore the collinear radiators
               IGNSOFT = 1    << 7,     // ignore the soft function
               IGNSNGL = 1    << 8,     // ignore the soft non-global function
               IGNPDF = 1     << 9,     // ignore the pdf ratio
               IGNFFUNC = 1   << 10,     // ignore ffunction
               EXPAND = COLLEXPAND | SOFTEXPAND | PDFEXPAND, // conveninece def.
               HYPGEO = 1     << 11,    // include hypergeometric function for transition point
    };

    const std::map<std::string,MODE> m_ModeToEnum = {{"RESUM", MODE::RESUM},
                                                     {"PDFEXPAND", MODE::PDFEXPAND},
                                                     {"COLLEXPAND", MODE::COLLEXPAND},
                                                     {"SOFTEXPAND", MODE::SOFTEXPAND},
                                                     {"LARGENC", MODE::LARGENC},
                                                     {"EXPAND", MODE::EXPAND},
                                                     {"HYPGEO", MODE::HYPGEO}, 
                                                     // from here on, modes are mostly usefull for debugging,
                                                     // should not appear in physics calculations
                                                     {"CKINV", MODE::CKINV},
                                                     {"CKCOUL", MODE::CKCOUL},
                                                     {"CHECKINVERSE", MODE::CKINV},
                                                     {"CHECKCOULOMB", MODE::CKCOUL},
                                                     {"CHECKCOULOMB", MODE::CKCOUL},
                                                     {"IGNCOLL", MODE::IGNCOLL},
                                                     {"IGNORECOLLINEAR", MODE::IGNCOLL},
                                                     {"IGNSOFT", MODE::IGNSOFT},
                                                     {"IGNORESOFT", MODE::IGNSOFT},
                                                     {"IGNSNGL", MODE::IGNSNGL},
                                                     {"IGNORESOFTNONGLOBAL", MODE::IGNSNGL},
                                                     {"IGNPDF", MODE::IGNPDF},
                                                     {"IGNOREPDF",MODE::IGNPDF},
                                                     {"IGNFFUNC",MODE::IGNFFUNC},
                                                     {"IGNOREFFUNCTION",MODE::IGNFFUNC}};

    enum MATCH_MODE {
      NONE = 0,
      LO   = 1       << 1,
      NLO = 1        << 2,
      ADD  = 1       << 3,
      DERIV = 1      << 4
    };

    const std::map<std::string,MATCH_MODE> m_MModeToEnum = {{"NONE", MATCH_MODE::NONE},
                                                            {"LO", MATCH_MODE::LO},
                                                            {"NLO", MATCH_MODE::NLO},
                                                            {"ADD", MATCH_MODE::ADD},
                                                            {"DERIV",MATCH_MODE::DERIV}};

      
    MODE m_amode = MODE::RESUM;
    size_t m_LogOrd = 1;

    // prinicipal groom mode
    RESUM::GROOM_MODE m_gmode = GROOM_MODE::NONE;
    // groom modes for legs, only used if m_gmode not NONE
    std::vector<RESUM::GROOM_MODE> m_collgmodes;
    RESUM::GROOM_MODE m_softgmode = GROOM_MODE::SD;

    std::valarray<std::vector<RESUM::GROOM_MODE>> m_allCollgmodes;
    std::valarray<RESUM::GROOM_MODE> m_allSoftgmodes;


    std::vector<RESUM::GROOM_MODE> m_collgmodes_end;
    RESUM::GROOM_MODE m_softgmode_end = GROOM_MODE::SD;

    MATCH_MODE m_mmode = MATCH_MODE::NONE;

    
    double CalcS(const double L, const double LResum, double &SoftexpNLL_LO, double& SoftexpNLL_NLO, MODE Check=MODE::RESUM);

    std::valarray<double> CalcS(std::valarray<double>& SoftexpNLL_LO, 
                                std::valarray<double>& SoftexpNLL_NLO, 
                                MODE Check=MODE::RESUM);

    double CalcColl(const double L, const double LResum, const int order,
                    double &Rp, MatrixD& G, MatrixD& Rexp, double& S1) {
      double dummy1;
      return CalcColl(L, LResum, order, Rp, G, Rexp, S1, dummy1);
    }
    double CalcColl_end(const double L, const double LResum, const int order,
                    double &Rp, MatrixD& G, MatrixD& Rexp, double& RAtEnd) {
      double dummy1;
      return CalcColl(L, LResum, order, Rp, G, Rexp, dummy1, RAtEnd);
    }
    double CalcColl(const double L, const double LResum, const int order,
                    double &Rp, MatrixD& G, MatrixD& Rexp, double& S1, double& RAtEnd);
    
    double CalcRpp(const double L, GROOM_MODE gmode, double &exp12);

    double CalcPDF(const double L, const double LResum, double &PDFexp);

    double CollinearCounterTerms( const int i,
                                  const ATOOLS::Flavour &fl,
                                  const ATOOLS::Vec4D &p,
                                  const double &z,
                                  const double muF2) const;

    double T(const double& x);
    std::valarray<double> T(const std::valarray<double>& x);

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

    virtual std::string AddObservable(const RESUM::Observable_Key& key,
                                      const std::valarray<double>& xvals);
   
    size_t nLegs() const {return n_g+n_q+n_aq+color_sings;}
    size_t nColoredLegs() const {return n_g+n_q+n_aq;}
    std::pair<int, double> Result(int n) const { return m_ress.at(n); }
    bool Initialized() const { return m_init; } 
    virtual void ResetObservables();

    std::vector<std::valarray<double>> m_resNLL;
    std::vector<double> resNLL(size_t i) { return {std::begin(m_resNLL[i]),std::end(m_resNLL[i])};}
    std::vector<std::valarray<double>> m_resExpLO;
    std::vector<double> resExpLO(size_t i) { return {std::begin(m_resExpLO[i]),std::end(m_resExpLO[i])};}
    std::vector<std::valarray<double>> m_resExpNLO;
    std::vector<double> resExpNLO(size_t i) { return {std::begin(m_resExpNLO[i]),std::end(m_resExpNLO[i])};}

    std::vector<std::valarray<double>> m_xvals;
    
    size_t m_nx;
    double m_ep;
    double m_p;
    double m_logFac;

    std::valarray<double> m_logArg;
    std::valarray<double> m_L;
    double m_Lz;
    std::valarray<double> m_Lsoft;
    std::valarray<double> m_epRatio;
    std::valarray<double> m_lambda;
    std::valarray<double> m_TofLoverA;


    double m_alphaS;
    double m_beta0;
    




  };// end of class Resum
  
}// end of namespace RESUM

#endif
