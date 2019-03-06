#include "AddOns/Analysis/Analyses/Analysis_Base.H"
#include "Analysis/Observable_Base.H"

namespace RESUM {

  class Matching_Analysis: public ANALYSIS::Analysis_Base {
  public:

    Matching_Analysis(const ANALYSIS::Argument_Matrix &params);

    ~Matching_Analysis();

    
    void Evaluate(double weight, double ncount,int mode);
    Primitive_Observable_Base * Copy() const;

    enum MODE {
               DEFAULT = 0,
               SKIP_REAL = 1        << 0,  // skip the real correction
               SKIP_SUBT = 1        << 1,  // skip all subtraction terms
               SKIP_SCOL = 1        << 2,  // skip the term related to the soft part of the splitting function
               SKIP_COLL = 1        << 3,  // skip the term related to the collinear part if the splitting function
               SKIP_NORM = 1        << 4,  // skip the term related to the normalization of the observable (dbar - b log(...))
               SKIP_SOFT = 1        << 5,  // skip the soft anomalous dimenstion matrix
               SKIP_PDFR = 1        << 6   // skip the pdf ratio
    };

    const std::map<std::string,MODE> m_ModeToEnum = {{"DEFAULT", MODE::DEFAULT},
                                                     {"SKIP_REAL", MODE::SKIP_REAL},
                                                     {"SKIP_SUBT", MODE::SKIP_SUBT},
                                                     {"SKIP_SCOL", MODE::SKIP_SCOL},
                                                     {"SKIP_COLL", MODE::SKIP_COLL},
                                                     {"SKIP_NORM", MODE::SKIP_NORM},
                                                     {"SKIP_SOFT", MODE::SKIP_SOFT},
                                                     {"SKIP_PDFR", MODE::SKIP_PDFR}};


    
  private:

    ANALYSIS::Argument_Matrix m_params;

    std::vector<Observable_Base*> m_obss;
    
    MODE m_mode;

    double EtaBeam(const ATOOLS::Vec4D& p, size_t beamId);
    

  };// end of class Matching_Analysis

}// end of namespace RESUM

#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "Tools/StringTools.H"
#include <algorithm>

#ifdef USING__ROOT
#undef USING__ROOT
#include "TH2D.h"
#include "AddOns/Analysis/Tools/My_Root.H"
#endif

using namespace RESUM;
using namespace ANALYSIS;
using namespace ATOOLS;

Matching_Analysis::Matching_Analysis(const Argument_Matrix &params):
  Analysis_Base(params[0][0]), m_params(params)
{
  DEBUG_FUNC(this);
  m_name+="_Resum";
  Data_Reader reader(",",";","!","=");
  const std::string& mode = reader.GetValue<std::string>("MATCHING_TEST","DEFAULT");
  if(is_int(mode)) {
    m_mode = static_cast<MODE>(to_type<int>(mode));
  }
  else {
    for(const std::string& m: split(mode,"\\|")) {
      m_mode = static_cast<MODE>(m_mode | m_ModeToEnum.at(m));
    }
  }
  Algebra_Interpreter *ip=reader.Interpreter();
  for (size_t i(1);i<params.size();++i) {
    if (params[i].size()<5) continue;
    Observable_Base *obs(Observable_Getter::GetObject
			 (params[i][0],Observable_Key(params[i][0])));
    if (obs == nullptr) {
      msg_Error()<<METHOD<<"(): Observable not found '"<<params[i][0]<<"'.\n";
      continue;
    }
    double xmin = ToType<double>(ip->Interprete(params[i][1]));
    double xmax = ToType<double>(ip->Interprete(params[i][2]));
    size_t nbin = ToType<size_t>(ip->Interprete(params[i][3]));
    int tp = HistogramType(params[i][4]);
    msg_Debugging()<<"Init '"<<params[i][0]<<"', type "<<tp
		   <<" with "<<nbin<<" bins in ["<<xmin<<","<<xmax<<"]\n";
    m_histos.push_back(new Histogram(tp,xmin,xmax,nbin,params[i][0]));
    m_obss.push_back(obs);
  }
  
#ifdef USING__ROOT
  if (MYROOT::myroot==NULL) {
    MYROOT::myroot = new MYROOT::My_Root();
    MYROOT::myroot->InitFile();
    ATOOLS::exh->AddTerminatorObject(MYROOT::myroot);
  }
  (*MYROOT::myroot)(new TH2D("yzplot","yzplot",100,-16.,0.,100,-16.,0.),"yzplot");
#endif

}

Matching_Analysis::~Matching_Analysis()
{
#ifdef USING__ROOT
  if (MYROOT::myroot) {
    exh->RemoveTerminatorObject(MYROOT::myroot);
    delete MYROOT::myroot;
    MYROOT::myroot=NULL;
  }
#endif
  for (size_t i(0);i<m_obss.size();++i) delete m_obss[i];
}

double Matching_Analysis::EtaBeam(const Vec4D& p, size_t beamId) {
  double eta = 1./(2.0*rpa->gen.PBeam(beamId)[0]);
  if(beamId == 0) {
    eta *= p.PPlus();
  }
  else if(beamId == 1) {
    eta *= p.PMinus();
  }
  else THROW(fatal_error, "Something went wrong, beamId > 1.");
  return eta;
}

void Matching_Analysis::Evaluate(double weight,double ncount,int mode)
{
  DEBUG_FUNC("mode = "<<mode);
  if (mode==0) {
    // calculate observable for event and fill into histo
    Particle_List all(*p_ana->GetParticleList(m_listname));
    Vec4D_Vector mom(2+all.size());
    Flavour_Vector fl(2+all.size());
    for (size_t i=0; i<all.size(); i++) {
      mom[2+i]=all[i]->Momentum();
      fl[2+i]=all[i]->Flav();
    }
    fl[0]=rpa->gen.Beam1();
    fl[1]=rpa->gen.Beam2();
    for (size_t i=0; i<m_obss.size(); i++) {
      double value = m_obss[i]->Value(&mom[0], &fl[0], mom.size());
      msg_Debugging()<<"value["<<i<<"] = "<<value
		     <<" ( w = "<<weight<<", n = "<<ncount<<" )\n";
      FillHisto(i,value,weight,ncount,mode);
    }
    return;
  }
  NLO_subevt* sub = p_ana->Sub();
  NLO_subevt* real = p_ana->Real();
  if (sub==real) {
    for (size_t i(0);i<m_obss.size();++i) {
      double value = m_obss[i]->Value(real->p_mom,real->p_fl,real->m_n);
      if (!(m_mode & MODE::SKIP_REAL)) FillHisto(i,value,weight,ncount,mode);
    }
    return;
  }
  if (m_mode & MODE::SKIP_SUBT) return;
  if (sub==nullptr) THROW(fatal_error,"Missing subtraction info.");
  for (size_t i=0; i<m_obss.size(); i++) {
    RESUM::Obs_Params ps = m_obss[i]->Parameters(sub->p_mom,sub->p_fl,
                                                 sub->m_n,sub->m_ijt);
    const double y = sub->m_y;
    double z = sub->m_z;
    double tau = y;

    double lrat = 0;
    if(!(m_mode & MODE::SKIP_COLL)) {
      lrat += sub->m_lt[1]/sub->m_lt[0];
    }
    if ((sub->m_i<2 && sub->p_fl[sub->m_ijt]==real->p_fl[sub->m_i]) ||
	(sub->m_i>=2 && (sub->p_fl[sub->m_ijt]==real->p_fl[sub->m_i] ||
			 sub->p_fl[sub->m_ijt]==real->p_fl[sub->m_j]))) {
      if (!(m_mode & MODE::SKIP_NORM)) {
        lrat+=2.0*m_obss[i]->Shift(sub)/sub->m_lt[0]/ps.m_a;
      }
      if (!(m_mode & MODE::SKIP_SOFT)) {
	double sij=dabs(2.0*sub->p_mom[sub->m_ijt]*sub->p_mom[sub->m_kt]);
        lrat+=log(sij/sub->m_mu2[stp::res])/sub->m_lt[0]*(ps.m_a+ps.m_b)/ps.m_a;
      }
    }

    // Calculate pdf ratio
    double pdfrat = 1.0;
    if (sub->m_i<2 && real->p_fl[sub->m_i].Strong()) {
      // Initial + something
      const double eta = abs(EtaBeam(sub->p_mom[sub->m_ijt], sub->m_i));
      z= (z-eta)/(1.0-eta);      

      if(!(m_mode & MODE::SKIP_PDFR)) {
        PHASIC::Single_Process *proc=static_cast<PHASIC::Single_Process*>
          ((void*)ToType<long unsigned int>(rpa->gen.Variable("CURPROC")));
        PDF::PDF_Base* pdf(proc->Integrator()->ISR()->PDF(sub->m_i));
        pdf->Calculate(eta, sub->m_mu2[stp::fac]);
        pdfrat *= pdf->GetXPDF(sub->p_fl[sub->m_ijt]);
        pdf->Calculate(eta/sub->m_z,sub->m_mu2[stp::fac]);
        pdfrat /= pdf->GetXPDF(real->p_fl[sub->m_i]);
      }

      if (sub->m_k<2) {
        // Initial + Initial
        tau /= 1.0-eta;
      }
      // else {
      // // Initial + Final
      // }
    }
    else if (sub->m_k<2 && real->p_fl[sub->m_k].Strong()) {
      // Final + Initial
      const double eta = abs(EtaBeam(sub->p_mom[sub->m_kt], sub->m_k));
      tau /= 1.0-eta;
      
      if(!(m_mode & MODE::SKIP_PDFR)) {
        PHASIC::Single_Process *proc=static_cast<PHASIC::Single_Process*>
          ((void*)ToType<long unsigned int>(rpa->gen.Variable("CURPROC")));
        PDF::PDF_Base* pdf(proc->Integrator()->ISR()->PDF(sub->m_k));
        pdf->Calculate(eta,sub->m_mu2[stp::fac]);
        pdfrat*=pdf->GetXPDF(sub->p_fl[sub->m_kt]);
        pdf->Calculate(eta/(1.0-sub->m_y),sub->m_mu2[stp::fac]);
        pdfrat/=pdf->GetXPDF(real->p_fl[sub->m_k]);
      }
    }
    
    if (!(m_mode & MODE::SKIP_SCOL)) {
      if (pow(tau,1.0/ps.m_a) < (1.0-z)) lrat+=sub->m_lt[2]/sub->m_lt[0];
      if (pow(tau,1.0/ps.m_a) < z) lrat+=sub->m_lt[3]/sub->m_lt[0];
    }
    lrat *= 2.0/(ps.m_a+ps.m_b);
    
#ifdef USING__ROOT
    ((TH2D*)(*MYROOT::myroot)["yzplot"])->Fill
      (log(tau),log(1.0-z),dabs(weight)*lrat*pdfrat);
#endif
    FillHisto(i, tau, weight*lrat*pdfrat, ncount, mode);
  }
}

Primitive_Observable_Base *Matching_Analysis::Copy() const 
{
  return new Matching_Analysis(m_params);
}

DECLARE_GETTER(Matching_Analysis,"MatchResum",
	       Primitive_Observable_Base,Argument_Matrix);

Primitive_Observable_Base *ATOOLS::Getter
<Primitive_Observable_Base,Argument_Matrix,Matching_Analysis>::
operator()(const Argument_Matrix &parameters) const
{
  if (parameters.size()==0 || parameters[0].size()<1) return NULL;
  return new Matching_Analysis(parameters);
}

void ATOOLS::Getter
<Primitive_Observable_Base,Argument_Matrix,Matching_Analysis>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"list"; 
}
