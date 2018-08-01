#ifndef RESUM_FFUNCTION_H
#define RESUM_FFUNCTION_H

#include "PDF/Main/Shower_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "Main/Comix_Interface.H"
#include "Main/Cluster_Definitions.H"
#include "Analysis/Observable_Base.H"
#include "Main/Resum.H"

#include "FFunction/MultiplePrecision.H"
#include "FFunction/Observables.H"

using std::vector;
namespace RESUM {

  class FFunction: public Resum {
  private:

    Comix_Interface *p_comix;

    Cluster_Definitions *p_clus;

    ATOOLS::Cluster_Amplitude *p_ampl;

    PDF::PDF_Base    ** p_pdf;
    
    MODEL::Running_AlphaS *p_as;
    
    //list of available and initialized bases
    CMetric_Map    m_cmetrics;
    //suitable basis for actual event
    CMetric_Base * p_cmetric;
    
    std::vector<Observable_Base_MP::Ptr>     m_obss;
    std::vector<ATOOLS::Histogram*>   m_hist;
    std::vector<std::vector<double> > m_ress;

    std::vector<ATOOLS::Vec4D> m_cms, momlabels;
    std::vector<ATOOLS::Flavour> flavlabels;
    std::vector<double> m_Qij, m_a, m_b, m_logdbar;
    std::vector<int> signlabels;
    int n_g, n_q, n_aq, color_sings, m_amode;
    double s_12, m_rn[2];


    double Value(const double &v);

    void MakeKinematics(const mpfloat& z, const mpfloat& y, const mpfloat& phi,
                        unsigned int ij, unsigned int k, vector<Vec4MP>& event);

    void GenerateEvent(double v, double Q2,
                       double epsilon, const mpfloat& rho,
                       vector<Vec4MP>& event);
  public:

    FFunction(PDF::ISR_Handler *const isr,
	  MODEL::Model_Base *const model);

    ~FFunction();

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

    size_t AddObservable(Observable_Base_MP::Ptr const obs,
			 ATOOLS::Histogram *const h);

    const std::vector<double> &Result(int n) const { return m_ress[n]; }


  };// end of class FFunction

}// end of namespace RESUM

#endif
