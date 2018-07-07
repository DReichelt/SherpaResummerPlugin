#include "PDF/Main/Shower_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "Main/Comix_Interface.H"
#include "Math/Matrix.H"
#include "Tools/CMetric_Base.H"


namespace RESUM {

  class RRatios: public PDF::Shower_Base {
  private:

    Comix_Interface m_comix;


    ATOOLS::Cluster_Amplitude *p_ampl_np1;
    ATOOLS::Cluster_Amplitude *p_ampl_n;
    // TODO: have this to keep some unused functions happy
    // ATOOLS::Cluster_Amplitude *p_ampl;

    
    ATOOLS::Cluster_Leg* p_emit = nullptr;
    ATOOLS::Cluster_Leg* p_spect = nullptr;
    ATOOLS::Cluster_Leg* p_soft = nullptr;
    ATOOLS::Cluster_Leg* p_emit_n = nullptr;
    ATOOLS::Cluster_Leg* p_spect_n = nullptr;

    std::vector<size_t> m_ordered_ids;
        
    
    MODEL::Running_AlphaS *p_as;
    
    //list of available and initialized bases
    CMetric_Map    m_cmetrics;
    //suitable basis for actual event
    CMetric_Base * p_cmetric_n;
    CMetric_Base * p_cmetric_np1;
    
    std::vector<ATOOLS::Histogram*>   m_hist;
    std::vector<std::vector<double> > m_ress;

    int n_g, n_q, n_aq, color_sings, m_amode;


  public:

    RRatios(PDF::ISR_Handler *const isr,
	  MODEL::Model_Base *const model);

    ~RRatios();

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

    const std::vector<double> &Result(int n) const { return m_ress[n]; }


  };// end of class RRatios

  inline void SetMomenta(ATOOLS::Cluster_Amplitude* ampl, const std::vector<ATOOLS::Vec4D>& moms) {
    // TODO: implement some range checks to get usefull errors
    for(size_t i=0; i<moms.size(); i++) {
      if(i<ampl->NIn()) ampl->Leg(i)->SetMom(-moms[i]);
      else ampl->Leg(i)->SetMom(moms[i]);
    }
  }

  
  
}// end of namespace RESUM
