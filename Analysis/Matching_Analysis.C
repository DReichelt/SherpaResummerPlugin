#include "AddOns/Analysis/Analyses/Analysis_Base.H"
#include "Analysis/Observable_Base.H"

namespace RESUM {

  class Matching_Analysis: public ANALYSIS::Analysis_Base {
  private:

    ANALYSIS::Argument_Matrix m_params;

    std::vector<Observable_Base*> m_obss;

    int m_mode;

  public:

    Matching_Analysis(const ANALYSIS::Argument_Matrix &params);

    ~Matching_Analysis();

    void Evaluate(double weight, double ncount,int mode);
    Primitive_Observable_Base * Copy() const;

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
  m_mode=reader.GetValue<int>("MATCHING_TEST",0);
  Algebra_Interpreter *ip=reader.Interpreter();
  for (size_t i(1);i<params.size();++i) {
    if (params[i].size()<5) continue;
    Observable_Base *obs(Observable_Getter::GetObject
			 (params[i][0],Observable_Key(params[i][0])));
    if (obs==NULL) {
      msg_Error()<<METHOD<<"(): Observable not found '"<<params[i][0]<<"'.\n";
      continue;
    }
    double xmin(ToType<double>(ip->Interprete(params[i][1])));
    double xmax(ToType<double>(ip->Interprete(params[i][2])));
    size_t nbin(ToType<size_t>(ip->Interprete(params[i][3])));
    int tp(HistogramType(params[i][4]));
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

void Matching_Analysis::Evaluate(double weight,double ncount,int mode)
{
  DEBUG_FUNC("mode = "<<mode);
  if (mode==0) {
    Particle_List all(*p_ana->GetParticleList(m_listname));
    Vec4D_Vector mom(2+all.size());
    Flavour_Vector fl(2+all.size());
    for (size_t i(0);i<all.size();++i) {
      mom[2+i]=all[i]->Momentum();
      fl[2+i]=all[i]->Flav();
    }
    fl[0]=rpa->gen.Beam1();
    fl[1]=rpa->gen.Beam2();
    for (size_t i(0);i<m_obss.size();++i) {
      double value(m_obss[i]->Value(&mom[0],&fl[0],mom.size()));
      msg_Debugging()<<"value["<<i<<"] = "<<value
		     <<" ( w = "<<weight<<", n = "<<ncount<<" )\n";
      FillHisto(i,value,weight,ncount,mode);
    }
    return;
  }
  NLO_subevt *sub(p_ana->Sub()), *real(p_ana->Real());
  msg_Debugging()<<*sub<<"\n";
  if (sub==real) {
    for (size_t i(0);i<m_obss.size();++i) {
      double value(m_obss[i]->Value(real->p_mom,real->p_fl,real->m_n));
      msg_Debugging()<<"value["<<i<<"] = "<<value<<"\n";
      if (!(m_mode&1)) FillHisto(i,value,weight,ncount,mode);
    }
    return;
  }
  if (m_mode&2) return;
  if (sub==NULL) THROW(fatal_error,"Missing subtraction info");
  for (size_t i(0);i<m_obss.size();++i) {
    RESUM::Obs_Params ps(m_obss[i]->Parameters(sub->p_mom,sub->p_fl,sub->m_n,sub->m_ijt));
    double y(sub->m_y), z(sub->m_z), tau(y);
    double lrat((m_mode&8)?0.0:sub->m_lt[1]/sub->m_lt[0]), pdfrat(1.0);
    if ((sub->m_i<2 && sub->p_fl[sub->m_ijt]==real->p_fl[sub->m_i]) ||
	(sub->m_i>=2 && (sub->p_fl[sub->m_ijt]==real->p_fl[sub->m_i] ||
			 sub->p_fl[sub->m_ijt]==real->p_fl[sub->m_j]))) {
      if (!(m_mode&16)) lrat+=2.0*m_obss[i]->Shift(sub)/sub->m_lt[0]/ps.m_a;
      if (!(m_mode&32)) {
	double sij=dabs(2.0*sub->p_mom[sub->m_ijt]*sub->p_mom[sub->m_kt]);
	lrat+=log(sij/sub->m_mu2[stp::res])/sub->m_lt[0]*(ps.m_a+ps.m_b)/ps.m_a;
      }
    }
    if (sub->m_i<2 && real->p_fl[sub->m_i].Strong()) {
      PHASIC::Single_Process *proc=static_cast<PHASIC::Single_Process*>
	((void*)ToType<long unsigned int>(rpa->gen.Variable("CURPROC")));
      PDF::PDF_Base *pdf(proc->Integrator()->ISR()->PDF(sub->m_i));
      Vec4D pa(sub->p_mom[sub->m_ijt]);
      double eta=dabs(sub->m_i?pa.PMinus():pa.PPlus())/(2.0*rpa->gen.PBeam(sub->m_i)[0]);
      pdf->Calculate(eta,sub->m_mu2[stp::fac]);
      pdfrat*=pdf->GetXPDF(sub->p_fl[sub->m_ijt]);
      pdf->Calculate(eta/sub->m_z,sub->m_mu2[stp::fac]);
      pdfrat/=pdf->GetXPDF(real->p_fl[sub->m_i]);
      z=(z-eta)/(1.0-eta);
      if (sub->m_k<2) tau/=1.0-eta;
    }
    else if (sub->m_k<2 && real->p_fl[sub->m_k].Strong()) {
      PHASIC::Single_Process *proc=static_cast<PHASIC::Single_Process*>
	((void*)ToType<long unsigned int>(rpa->gen.Variable("CURPROC")));
      PDF::PDF_Base *pdf(proc->Integrator()->ISR()->PDF(sub->m_k));
      Vec4D pa(sub->p_mom[sub->m_kt]);
      double eta=dabs(sub->m_k?pa.PMinus():pa.PPlus())/(2.0*rpa->gen.PBeam(sub->m_k)[0]);
      pdf->Calculate(eta,sub->m_mu2[stp::fac]);
      pdfrat*=pdf->GetXPDF(sub->p_fl[sub->m_kt]);
      pdf->Calculate(eta/(1.0-sub->m_y),sub->m_mu2[stp::fac]);
      pdfrat/=pdf->GetXPDF(real->p_fl[sub->m_k]);
      tau/=1.0-eta;
    }
    if (!(m_mode&4)) {
    if (pow(tau,1.0/ps.m_a)<(1.0-z)) lrat+=sub->m_lt[2]/sub->m_lt[0];
    if (pow(tau,1.0/ps.m_a)<z) lrat+=sub->m_lt[3]/sub->m_lt[0];
    }
    lrat*=2.0/(ps.m_a+ps.m_b);
#ifdef USING__ROOT
    ((TH2D*)(*MYROOT::myroot)["yzplot"])->Fill
      (log(tau),log(1.0-z),dabs(weight)*lrat*pdfrat);
#endif
    FillHisto(i,tau,weight*lrat*pdfrat,ncount,mode);
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
