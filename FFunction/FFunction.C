#include "FFunction/FFunction.H"

#include "PHASIC++/Process/Single_Process.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Matrix.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Histogram.H"
#include "Analysis/Observable_Base.H"
#include "Math/matexp.hpp"



#include <typeinfo>

using namespace RESUM;
using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;
using namespace MODEL;

FFunction::FFunction(ISR_Handler *const isr,
	     Model_Base *const model):
Resum(isr,model), p_ampl(NULL)
{
  p_comix = new Comix_Interface();
  p_clus = new Cluster_Definitions();
  p_as=(Running_AlphaS*)model->GetScalarFunction("alpha_S");
  
  p_pdf = new PDF_Base*[2];
  for (int i=0;i<2; i++) p_pdf[i] = isr->PDF(i);

  Data_Reader read(" ",";","#","=");
  m_amode=read.GetValue<int>("RESUM_MODE",0);
  if (rpa->gen.Variable("SHOWER_GENERATOR")=="")
    rpa->gen.SetVariable("SHOWER_GENERATOR",ToString(this));
}

FFunction::~FFunction()
{
  if (p_ampl) p_ampl->Delete();
  delete p_comix;
  delete p_clus;
  delete [] p_pdf;

  while (m_cmetrics.size()>0) {
    delete m_cmetrics.begin()->second;
    m_cmetrics.erase(m_cmetrics.begin());
  }
}


Vec4MP Cross(const Vec4MP& a, const Vec4MP& b) {
  return {0, a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[3]};
}

Vec4MP Boost(const Vec4MP& a, const Vec4MP& b) {
  mpfloat rsq = a.Abs();
  mpfloat v0 = a*b/rsq;
  mpfloat c1 = (b[0]+v0)/(rsq+a[0]);
  return {v0, b[1]-c1*a[1], b[2]-c1*a[2], b[3]-c1*a[3]};
}

Vec4MP BoostBack(const Vec4MP& a, const Vec4MP& b) {
  mpfloat rsq = a.Abs();
  mpfloat v0 = (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]*a[3]*b[3])/rsq;
  mpfloat c1 = (b[0]+v0)/(rsq+a[0]);
  return {v0, b[1]+c1*a[1], b[2]+c1*a[2], b[3]+c1*a[3]};
  
}


void FFunction::MakeKinematics(const mpfloat& z, const mpfloat& y, const mpfloat& phi,
                    unsigned int ij, unsigned int k, vector<Vec4MP>& event){ 
  const auto& pijt = event.at(ij);
  const auto& pkt = event.at(k);
  Vec4MP Q = pijt+pkt;
  mpfloat rkt = mp::sqrt(Q.Abs2()*y*z*(1.-z));
  Vec4MP kt1 = Cross(pijt, pkt);
  if(kt1.P() < 1.e-6) kt1 =  Cross(pijt, Vec4MP(0.,1.,0.,0.));
  kt1 *= rkt*mp::cos(phi)/kt1.P();
  Vec4MP kt2cms =  Cross(Boost(Q,pijt), kt1);
  kt2cms *= rkt*mp::sin(phi)/kt2cms.P();
  Vec4MP kt2 = BoostBack(Q, kt2cms);
  mpfloat oneMz = (1-z);
  mpfloat oneMzy =  oneMz*y;
  mpfloat zy = z*y;
  mpfloat oneMy = mpfloat(1)-y;
  Vec4MP pij = {pijt[0],pijt[1],pijt[2],pijt[3]};
  Vec4MP pk = {pkt[0],pkt[1],pkt[2],pkt[3]};
  event.at(ij) = z*pij + oneMzy*pk + kt1 + kt2;
  event.at(k) = oneMy*pk;
  event.emplace_back(oneMz*pij + zy*pk - kt1 - kt2);
}



void FFunction::GenerateEvent(double v, double Q2,
                              double epsilon, const mpfloat& rho,
                              vector<Vec4MP>& event) {
  
  if(v < 1) {
  
    double CF = 4./3.;
  
    double L = -std::log(v);
    double asMax = p_as->AlphaS(v*v*Q2);
    double rpMax = 4*CF*asMax/M_PI * L;
    mpfloat Q = std::sqrt(Q2);
  
    double t = v*Q2;
    while(t > epsilon*v*Q2) {
      t *= std::pow(ran->Get(),1./rpMax);
      if(t > epsilon*v*Q2 && t < v*Q2) {
        double z = 1-std::exp(-ran->Get()*L);
        double weight = 1./2. * p_as->AlphaS((1.-z)*v*Q2)/asMax;
        
        if(weight > ran->Get()) {
          double kT2 = (1-z)*t;
          double eta = std::log((1-z)*std::sqrt(Q2/kT2));
          double phi = 2*M_PI*ran->Get();
          
          double xi = eta/(std::log(Q2/t)/2);
          mpfloat kT_rho = std::sqrt(kT2)*mp::pow(rho, 1-xi/2.);
          mpfloat eta_rho = eta-xi*mp::log(rho)/2.;
          mpfloat z_rho = mpfloat(1)-kT_rho/Q * mp::exp(eta_rho);
          
          unsigned int ij = ran->Get() > 0.5 ? 0 : 1;
          unsigned int k = ij == 0 ? 1 : 0;
          
          
          MakeKinematics(z_rho,
                         mp::pow(kT_rho/Q,2)/(z_rho*(1-z_rho)),
                         phi, ij+2, k+2, event);
        }
      }
    }
  }
}



int FFunction::PerformShowers()
{
  DEBUG_FUNC(this);
  if (p_ampl==NULL) THROW(fatal_error,"No process info");
  double epsilon = 0.000001;
  mpfloat rho = mp::pow(mpfloat(10),-500);


  msg_Debugging()<<"Using epsilon = "<<epsilon<<" and rho = "<<rho<<"\n";

  msg_Debugging()<<"No. of registered observables: "<<m_obss.size()<<"\n";
  

  if(m_obss.size()>0) {
msg_Debugging()<<"Observable: "<<m_obss.at(0)->Name()<<"\n";
    m_weight=1.0;
  vector<Vec4MP> moms_born(p_ampl->Legs().size());
  vector<Flavour> flavs_born(p_ampl->Legs().size());
  for (size_t i(0);i<p_ampl->Legs().size();++i) {
    moms_born[i]=i<p_ampl->NIn()?-p_ampl->Leg(i)->Mom():p_ampl->Leg(i)->Mom();
    flavs_born[i]=i<p_ampl->NIn()?p_ampl->Leg(i)->Flav().Bar():p_ampl->Leg(i)->Flav();
  }

  //TODO: loop over observables
  vector<Obs_Params> params(moms.size());
  for(size_t l=0; l<params.size(); l++) {
    params[l] = m_obss.at(0)->Parameters(moms, flavs, l);
  }

  //TODO: why choose only one bin at random?
  /* size_t i=1+m_hist.at(0)->Nbin()*ran->Get(); */
  /* double v = m_hist.at(0)->HighEdge(i); */

  for(size_t i=0; i<m_hist.at(0)->Nbin(); i++) {
    vector<Vec4MP> moms = moms_born;
    vector<Flavour> flavs = flavs_born;
    GenerateEvent(v, 91.2*91.2, epsilon, rho, moms);
    while(flavs.size()<moms.size()) flavs.push_back(21);
    m_obss.at(0)->Value(moms,flavs,p_ampl->NIn());
  }
  CleanUp();
  return 1;
}


int FFunction::PerformDecayShowers()
{
  DEBUG_FUNC(this);
  return PerformShowers();
}

bool FFunction::ExtractPartons(Blob_List *const bl)
{
  Blob *b(bl->FindLast(btp::Shower));
  if (b==NULL) THROW(fatal_error,"No Shower blob");
  b->SetTypeSpec("RESUM");
  for (int i=0;i<b->NInP();++i) 
    b->InParticle(i)->SetStatus(part_status::decayed);
  for (int i=0;i<b->NOutP();++i) 
    b->OutParticle(i)->SetStatus(part_status::decayed);
    b->SetStatus(blob_status::needs_beams |
	       blob_status::needs_hadronization);
  return true;
}

void FFunction::CleanUp()
{
  p_comix->Reset();
  if (p_ampl) p_ampl->Delete();
  p_ampl=NULL;
  
  m_cms.clear();
  m_Qij.clear();
  flavlabels.clear();
  signlabels.clear();
  momlabels.clear();

  m_a.clear();
  m_b.clear();
  m_logdbar.clear();
}

Cluster_Definitions_Base *FFunction::GetClusterDefinitions()
{
  return NULL;
}

bool FFunction::PrepareShower
(Cluster_Amplitude* ampl,const bool &soft)
{
  DEBUG_FUNC(this);
  DEBUG_VAR(ampl->Proc<Process_Base>());
  p_ampl = ampl->Copy();

  
  return true;
} 
  
size_t FFunction::AddObservable(Observable_Base_MP::Ptr const obs,
			    ATOOLS::Histogram *const h)
{
  m_obss.push_back(obs);
  m_hist.push_back(h);
  m_ress.push_back(std::vector<double>(2));
  return m_ress.size()-1;
}

double FFunction::CplFac(const ATOOLS::Flavour &fli,const ATOOLS::Flavour &flj,
		    const ATOOLS::Flavour &flk,const int type,
		    const int cpl,const double &mu2) const
{
  THROW(not_implemented,"");
  return -1.0;
}



DECLARE_GETTER(FFunction,"FFunction",Shower_Base,Shower_Key);

Shower_Base *Getter<Shower_Base,Shower_Key,FFunction>::
operator()(const Shower_Key &key) const
{
  return new FFunction(key.p_isr,key.p_model);
}

void Getter<Shower_Base,Shower_Key,FFunction>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"The MC for the FFunction."; 
}
