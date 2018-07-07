#include "PDF/Main/Shower_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "Main/Comix_Interface.H"
#include "Main/Cluster_Definitions.H"
#include "Analysis/Observable_Base.H"
#include "Tools/CMetric_Base.H"
#include "Tools/CBasis.H"
#include "Tools/CBasis.C"

namespace RESUM {

  class Resum: public PDF::Shower_Base {
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
    
    std::vector<Observable_Base*>     m_obss;
    std::vector<ATOOLS::Histogram*>   m_hist;
    std::vector<std::vector<double> > m_ress;

    std::vector<ATOOLS::Vec4D> m_cms, momlabels;
    std::vector<ATOOLS::Flavour> flavlabels;
    std::vector<double> m_Qij, m_a, m_b, m_logdbar;
    std::vector<int> signlabels;
    int n_g, n_q, n_aq, color_sings, m_amode;
    double s_12, m_rn[2];

    double CalcS(const double L,double &Softexp);
    double CalcF(const double L);
    double CalcColl(const double L,const int order,
		    double &Rp,double &Collexp);
    double CalcPDF(const double L,double &PDFexp);
    double T(const double x);
    void printMatrix(const Complex rhs[],const size_t dim);
    void matMult(Complex ResMat[],const Complex ArrMat[],
		 const std::vector<std::vector<double> > VecMat);
    Complex* H_conjugate(const Complex ResMat[],size_t dim);

    double Value(const double &v);

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

    size_t AddObservable(Observable_Base *const obs,
			 ATOOLS::Histogram *const h);

    const std::vector<double> &Result(int n) const { return m_ress[n]; }


  };// end of class Resum

    void Resum::printMatrix(const Complex rhs[], const size_t dim){
    msg_Debugging()<<"{\n";
    for (size_t i=0;i<dim;i++) {
      msg_Debugging()<<"  {";
      for (size_t j=0;j<dim;j++) {
	 msg_Debugging()<< (fabs(rhs[i*dim+j].real()) > .001 ? rhs[i*dim+j].real() : 0) <<"+"
		  <<rhs[i*dim+j].imag()<<((j+1<dim)?"I,":"I");
      }
      msg_Debugging()<<((i+1==dim)?"}\n":"},\n");
    }
    msg_Debugging()<<"}\n";  
}


//Matrix Multiplication
void Resum::matMult(Complex ResMat[], const Complex ArrMat[], const std::vector< std::vector< double > > VecMat){
    size_t dim = VecMat.size();
    for (size_t i=0;i<dim;i++) {
        for (size_t j=0;j<dim;j++){
    	   for (size_t k=0;k<dim;k++){
	     ResMat[i*dim+j] = Complex(ResMat[i*dim+j].real() + VecMat[k][i]*ArrMat[j*dim+k].real(),
				       ResMat[i*dim+j].imag() + VecMat[k][i]*ArrMat[j*dim+k].imag());
	   }      
	}
    }
}

//Hermition conjugate 
Complex* Resum::H_conjugate(const Complex ResMat[], size_t dim) {
complex<double> hold;
Complex *result = new Complex[dim*dim];
for (size_t i=0;i<dim;i++) {
        for (size_t j=0;j<dim;j++){ 
	  hold = ResMat[i*dim+j];
	  result[j*dim+i] = conj(hold);
	}
    }
 return result;
}

  
}// end of namespace RESUM
