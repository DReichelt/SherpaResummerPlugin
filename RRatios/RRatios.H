#include "PDF/Main/Shower_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "Main/Comix_Interface.H"
#include "Math/Matrix.H"
#include "Tools/CMetric_Base.H"


namespace RESUM {


  /// Compute ratio M_n*eikonal/M_n+1
  /**
     Class to automatically compute the ratio of the approximate 
     and full n+1 matrix element. Used to validate input files.
     Assumes that the rund is started with the n+1 matrix element
     as process, which hence has to include at least one gluon in 
     the final state.
   */
  class RRatios: public PDF::Shower_Base {
  private:

    /// The interface to comix, to compute hard matrix.
    Comix_Interface m_comix;

    /// Pointer to amplitude with n+1 partons.
    ATOOLS::Cluster_Amplitude *p_ampl_np1;

    /// Pointer to amplitude with n partons.
    ATOOLS::Cluster_Amplitude *p_ampl_n;
    // TODO: have this to keep some unused functions happy
    // ATOOLS::Cluster_Amplitude *p_ampl;

    /// Pointer to emitting leg in n+1-parton amplitude.
    ATOOLS::Cluster_Leg* p_emit = nullptr;

    /// Pointer to spectating leg taking the recoil in n+1-parton amplitude.
    ATOOLS::Cluster_Leg* p_spect = nullptr;

    /// Pointer to soft leg in n+1 amplitude.
    ATOOLS::Cluster_Leg* p_soft = nullptr;

    /// Pointer to emitter in n-parton amplitude.
    ATOOLS::Cluster_Leg* p_emit_n = nullptr;

    /// Pointer to spectating leg taking the recoil in n-parton amplitude.
    ATOOLS::Cluster_Leg* p_spect_n = nullptr;

    /// Leg ids ordered as quarks-antiquarks-gluons in n-parton amplitude.
    std::vector<size_t> m_ordered_ids;

    /// Leg ids ordered as quarks-antiquarks-gluons in n+1-parton amplitude.
    std::vector<size_t> m_ordered_ids_np1;
        
    /// Pointer to alphaS function.
    MODEL::Running_AlphaS *p_as;
    
    /// List of available and initialized bases.
    CMetric_Map    m_cmetrics;
    
    /// Suitable basis for actual n-parton amplitude.
    CMetric_Base * p_cmetric_n;

    /// Suitable basis for actual n+1-parton amplitude.
    CMetric_Base * p_cmetric_np1;

    /// \todo We dont really need this here.
    std::vector<ATOOLS::Histogram*>   m_hist;
    std::vector<std::vector<double> > m_ress;

    /// Number of gluons, quarks, colorsigns.
    /// \todo Why are these members?
    int n_g, n_q, n_aq, color_sings;

    /// Resum mode
    /// \todo Do we need this for RRatios?
    int m_amode;

    /// Center of mass Energy (-1 if not set).
    double m_E = -1;

    /// Multiplicative parameter steering how fast soft limit is approached.
    double m_lambda = 0.95;

    /// Cutoff for lambda
    double m_cutoff = 1e-3;
    
    /// Number of events handled so far. Used for names of plots.
    size_t m_count=0;

    /// Do or do not write ratio plot to file.
    bool m_plotRatio = true;
    
    /// Do or do not write matrix element for n process to file.
    bool m_plotMEn = true;
    
    /// Do or do not write matrix element for n+1 process to file.
    bool m_plotMEnp1 = true;
    
    /// Set the final state to equally distanced partons in the x-y plane.
    /// \todo Do not hardcode energy
    void MinimallyCollinearFinalState(ATOOLS::Cluster_Amplitude* ampl);

  public:

    /// \name Constructors and Destructors.
    // \{
    /// Standard and only constructor.
    RRatios(PDF::ISR_Handler *const isr,
	  MODEL::Model_Base *const model);

    /// Destructor.
    ~RRatios();
    // \}

    /// Specialisation of virtual function.
    int  PerformShowers();

    /// Specialisation of virtual function
    int  PerformDecayShowers();

    /// Specialisation of virtual function
    bool ExtractPartons(ATOOLS::Blob_List *const bl);

    /// Reset everything and do some cleanup. Legacy from Resum.
    void CleanUp();

    /// Specialisation of virtual function.
    /**
       Required by vitual class but just returns nullptr.
     */
    PDF::Cluster_Definitions_Base *GetClusterDefinitions();

    /// Specialisation of virtual function.
    /**
       Prepare everything for the calculation. In this case in particular
       identify emitter, spectator and soft gluon, set up the n-parton
       amplitude and find bases for both of them.
     */
    bool PrepareShower(ATOOLS::Cluster_Amplitude *const ampl,
		       const bool & soft=false);

    /// Specialisation of virtual function.
    /**
       Required by framework by not really used here.
     */
    double CplFac(const ATOOLS::Flavour &fli,const ATOOLS::Flavour &flj,
		  const ATOOLS::Flavour &flk,const int type,
		  const int cpl,const double &mu2) const;

    /// Accesor to m_ress.
    /**
       \todo This is required by the analysis stuff that is not really used here anyway.
     */
    const std::vector<double> &Result(int n) const { return m_ress[n]; }


  };// end of class RRatios

  /// Set momenta properly in ATOOLS::Cluster_Amplitude
  /**
     \todo This could probably be moved to some common "deal with sherpa conventions" place...
   */
  inline void SetMomenta(ATOOLS::Cluster_Amplitude* ampl, const std::vector<ATOOLS::Vec4D>& moms) {
    /// \todo Implement some range checks to get usefull errors.
    for(size_t i=0; i<moms.size(); i++) {
      if(i<ampl->NIn()) ampl->Leg(i)->SetMom(-moms[i]);
      else ampl->Leg(i)->SetMom(moms[i]);
    }
  }
}// end of namespace RESUM
