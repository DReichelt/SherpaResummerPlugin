#ifndef RESUM_RESUM_H
#define RESUM_RESUM_H

#include "PDF/Main/Shower_Base.H"

#include "Main/Comix_Interface.H"
#include "Main/Params.H"
#include "Tools/CMetric_Base.H"
#include "Analysis/Observable_Base.H"

#include "Analysis/Cumulant.H"
#include "Analysis/ChannelAlgorithms/ChannelAlgorithm_Base.H"

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
    std::map<std::string,Params> m_varParams;


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

    std::vector<std::valarray<RESUM::GROOM_MODE>> m_allCollgmodes;
    std::valarray<RESUM::GROOM_MODE> m_allSoftgmodes;


    std::vector<RESUM::GROOM_MODE> m_collgmodes_end;
    RESUM::GROOM_MODE m_softgmode_end = GROOM_MODE::SD;

    std::vector<std::valarray<bool>> m_collGroomed;
    std::valarray<bool> m_softGroomed;

    MATCH_MODE m_mmode = MATCH_MODE::NONE;

    

    std::valarray<double> CalcS(double& SoftexpNLL_LO, 
                                double& SoftexpNLL_NLO, 
                                MODE Check=MODE::RESUM);
    std::valarray<double> CalcSNGL(double& SnglExpNLL_NLO);
    std::valarray<double> CalcF(double& FexpNLL_NLO);
    std::valarray<double> CalcEP(std::valarray<double>& EPexpNLL_LO, std::valarray<double>& EPexpNLL_NLO);
    std::valarray<double> CalcColl(std::valarray<double>& Rp, Matrix<std::valarray<double>>& G, 
                                   Matrix<std::valarray<double>>& Rexp, double& S1, 
                                   std::valarray<double>& RAtEnd);
    double CalcRpp(const double L, GROOM_MODE gmode, double &exp12);
    std::valarray<double> CalcPDF(double &PDFexp);

    double CollinearCounterTerms( const int i,
                                  const ATOOLS::Flavour &fl,
                                  const ATOOLS::Vec4D &p,
                                  const double &z,
                                  const double muF2) const;

    double T(const double& x) const;
    std::valarray<double> T(const std::valarray<double>& x) const;
    
  public:

    Resum(PDF::ISR_Handler *const isr,
	  MODEL::Model_Base *const model);

    ~Resum();


    int PerformShowers(double weight) {m_weight=weight; return PerformShowers();}
    int PerformShowers();
    int PerformDecayShowers();

    bool ExtractPartons(ATOOLS::Blob_List *const bl);
    void CleanUp();

    PDF::Cluster_Definitions_Base *GetClusterDefinitions();

    bool PrepareShower(ATOOLS::Cluster_Amplitude *const ampl,
		       const bool & soft=false);
   
    void FinishShowers(double weight, int ntrials);

    double CplFac(const ATOOLS::Flavour &fli,const ATOOLS::Flavour &flj,
		  const ATOOLS::Flavour &flk,const int type,
		  const int cpl,const double &mu2) const;

    virtual size_t AddObservable(Observable_Base *const obs,
                                 ATOOLS::Histogram *const h);

    virtual std::string AddObservable(const RESUM::Observable_Key& key,
                                      std::valarray<double>& xvals);

    virtual std::string AddObservable(const RESUM::Observable_Key& key,
                                      std::valarray<double>& xvals,
                                      std::string title, bool calcNLL, 
                                      bool calcLO, bool calcNLO);
    
    void AddChAlg(const std::string& name, const std::vector<std::string> args);
    std::vector<std::string> m_channels;
    std::vector<std::string> m_curChannels;
    std::vector<ChannelAlgorithm_Base::Ptr> m_channelAlgs;
   
    size_t nLegs() const {return n_g+n_q+n_aq+color_sings;}
    size_t nColoredLegs() const {return n_g+n_q+n_aq;}
    std::pair<int, double> Result(int n) const { return m_ress.at(n); }
    bool Initialized() const { return m_init; } 
    virtual void ResetObservables();

    std::vector<std::string> m_varNames;
    const std::vector<std::string>& VarNames() {return m_varNames;}
    int m_varIdx = -1;

    void SetVariationParameters(const std::vector<SHERPA::Variation_Parameters *> * pvec);
    void SetVariationParameters(const SHERPA::Variations* vars);
    void SetVariationIdx(int idx) {m_varIdx=idx;}

    size_t VarNum() {return p_variationweights ? p_variationweights->GetNumberOfVariations() : 0;}
    int VarSubNum() {return p_variationweights ? p_variationweights->GetNumberOfSubevents() : 0;} 
    int VarIdx() {return std::min(m_varIdx,VarSubNum());}
    std::string VarName(size_t i) {return p_variationweights ? p_variationweights->GetVariationNameAt(i) : "";}
    double VarWeight(size_t i) {return p_variationweights ? p_variationweights->GetVariationWeightAt(i, SHERPA::Variations_Type::all, VarIdx()) : 0;}

    std::map<std::string, double> VarWeightsMap() {
      std::map<std::string, double> varweights;
      for(size_t i=0; i<VarNum(); i++) {
        varweights[VarName(i)] = VarWeight(i);
      }
      return varweights;
    }

    std::vector<std::valarray<double>> m_resNLL;
    std::vector<std::map<std::string,std::valarray<double>>> m_varResNLL;
    std::vector<double> resNLL(size_t i) { return {std::begin(m_resNLL[i]),std::end(m_resNLL[i])};}
    std::vector<std::valarray<double>> m_resExpLO;
    std::vector<std::map<std::string,std::valarray<double>>> m_varResExpLO;
    std::vector<double> resExpLO(size_t i) { return {std::begin(m_resExpLO[i]),std::end(m_resExpLO[i])};}
    std::vector<std::valarray<double>> m_resExpNLO;
    std::vector<std::map<std::string,std::valarray<double>>> m_varResExpNLO;
    std::vector<double> resExpNLO(size_t i) { return {std::begin(m_resExpNLO[i]),std::end(m_resExpNLO[i])};}

    std::set<size_t> m_orders;
    void SetOrders(std::set<size_t> orders) {m_orders = orders;}
    void AddOrder(size_t order) {m_orders.insert(order);}
    bool HasOrder(size_t order) {return m_orders.find(order) != m_orders.end();}

    std::map<std::string, Cumulant::Ptr> m_resultsNLL;
    std::map<std::string, Cumulant::Ptr> m_resultsExpLO;
    std::map<std::string, Cumulant::Ptr> m_resultsExpNLO;
    std::map<std::string, Cumulant::Ptr> m_resultsLO;
    std::map<std::string, Cumulant::Ptr> m_resultsNLO;


    enum TREAT {
      none = 0,
      born   = 1      << 1,
      lo = 1          << 2,
      nlo  = 1        << 3,
      needs_resum = 1 << 4
    };    

    std::vector<size_t> m_treat;
    size_t m_globaltreat = TREAT::none;
    std::vector<double> m_resFO;

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
    std::valarray<double> m_lambdaSoft;
    std::valarray<double> m_TofLoverA;
    std::valarray<double> m_Rp;

    // std::valarray<MatrixD> m_G;
    Matrix<std::valarray<double>> m_G;
    Matrix<std::valarray<double>> m_Rexp;
    double m_S1; 
    double m_S2;
    double m_P1;
    double m_F2;
    double m_SNGL2;
    std::valarray<double> m_EP1;
    std::valarray<double> m_EP2;
    std::valarray<double> m_RAtEnd;
      
    std::valarray<double> m_Coll;
    std::valarray<double> m_Soft;
    std::valarray<double> m_PDF;
    std::valarray<double> m_Ffunc;
    std::valarray<double> m_SNGL;
    std::valarray<double> m_Ep;


    double NC() const {return m_params.NC();}
    double CA() const {return m_params.CA();}
    double CF() const {return m_params.CF();}
    double TR() const {return m_params.TR();}
    bool largeNC() const {return m_params.largeNC();}
    double beta0() const {return m_params.beta0(p_ampl->MuR2());}
    double beta1() const {return m_params.beta1(p_ampl->MuR2());}
    double K_CMW() const {return m_params.K_CMW(p_ampl->MuR2());}
    double CollDimGlue() const {return m_params.CollDimGlue(p_ampl->MuR2());}
    double CollDimQuark() const {return m_params.CollDimQuark(p_ampl->MuR2());}
    double alphaS() const {return m_params.alphaS(p_ampl->MuR2());}
    double alphaSBar() const {return alphaS()/2./M_PI;}
    double muR2() const {return m_params.MuR2fac()*p_ampl->MuR2();}
    double muF2() const {return m_params.MuF2fac()*p_ampl->MuF2();}
    double muQ2() const {return m_obss[m_n]->MuQ2fac()*p_ampl->MuQ2();}
    double muR() const {return sqrt(muR2());}
    double muF() const {return sqrt(muF2());}
    double muQ() const {return sqrt(muQ2());}

    std::vector<double> m_energy;

    inline bool isQuark(size_t l) const {return p_ampl->Leg(l)->Flav().IsQuark();}
    inline bool isGluon(size_t l) const {return p_ampl->Leg(l)->Flav().IsGluon();}
    inline bool isStrong(size_t l) const {return isQuark(l) or isGluon(l);}

    inline double colfac(size_t l) const {return isQuark(l) ? CF() : (isGluon(l)? CA():0);}
    inline double C(size_t l) const {return colfac(l);}
    inline double hardcoll(size_t l) const {return isQuark(l) ? CollDimQuark() : (isGluon(l)? CollDimGlue():0);}
    inline double B(size_t l) const {return hardcoll(l);}
    inline double energy(size_t l) const {return m_energy[l];}
    inline double E(size_t l) const {return energy(l);}

    std::valarray<double> r1(std::valarray<double> lambda, double a, double b) const;
    std::valarray<double> r1_b0(std::valarray<double> lambda, double a) const;
    std::valarray<double> SD_r1(std::valarray<double> lambda, double lambdaZ, double a, double b) const;
    std::valarray<double> SD_r1_b0(std::valarray<double> lambda, double lambdaZ, double a) const;
    std::valarray<double> r2_cmw(std::valarray<double> lambda, double a, double b) const;
    std::valarray<double> r2_beta1(std::valarray<double> lambda, double a, double b) const;
    std::valarray<double> r2_cmw_b0(std::valarray<double> lambda, double a) const;
    std::valarray<double> r2_beta1_b0(std::valarray<double> lambda, double a) const;
    std::valarray<double> r1p_b0(std::valarray<double> lambda, double a) const;
    std::valarray<double> r1p(std::valarray<double> lambda, double a, double b) const;
    std::valarray<double> SD_r2_cmw(std::valarray<double> lambda, double lambdaZ, double a, double b) const;
    std::valarray<double> SD_r2_beta1(std::valarray<double> lambda, double lambdaZ, double a, double b) const;
    std::valarray<double> SD_r1p(std::valarray<double> lambda, double lambdaZ, double a, double b) const;
    std::valarray<double> SD_r1d(std::valarray<double> lambda, double lambdaZ, double a, double b) const;
    std::valarray<double> SD_r2_cmw_b0(std::valarray<double> lambda, double lambdaZ, double a) const;
    std::valarray<double> SD_r2_beta1_b0(std::valarray<double> lambda, double lambdaZ, double a) const;
    std::valarray<double> SD_r1p_b0(std::valarray<double> lambda, double lambdaZ, double a) const;
    std::valarray<double> SD_r1d_b0(std::valarray<double> lambda, double lambdaZ, double a) const;

    void printLists(std::ostream& stream, std::vector<std::valarray<double>> lists) {
      for(size_t i=0; i<m_nx; i++) {
        for(auto& l: lists) stream<<l[i]<<" ";
        stream<<"\n";
      }
    }
    

  };// end of class Resum
  
  double NaNtoZero(double x) {return std::isnan(x)?0.:x;}
  double FuzzyZeroToZero(double x) {return ATOOLS::IsZero(x)?0.:x;}
    
}// end of namespace RESUM

#endif
