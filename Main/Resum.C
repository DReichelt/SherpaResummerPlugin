#include "Main/Resum.H"

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



const double s_Nc = 3.;
const double s_CF = (s_Nc*s_Nc-1.)/(2.*s_Nc);
const double s_CA = s_Nc;
const double s_TR = 1./2.;
const double s_eps = .000001;


//Print a matrix
template <typename T>
void printMat(std::vector< std::vector <T> > &rhs) {
  msg_Debugging() << std::endl;
  for (unsigned i=0; i<rhs.size(); i++) {
    msg_Debugging() << "{";
    for (unsigned j=0; j<rhs[i].size(); j++) {
      msg_Debugging() << (fabs(rhs[i][j]) > .0000001 ? rhs[i][j] : 0);
      if(j<rhs[j].size()-1) msg_Debugging() << ",";
    }
    msg_Debugging() << "},";
    msg_Debugging() << std::endl;
  }
  msg_Debugging() << std::endl;
}
 
 
// returns minus a matrix
inline void Minus(std::vector< std::vector <double> > &ts){
  for(unsigned i = 0; i<ts.size(); i++){
    for(unsigned j = 0; j<ts[i].size(); j++){
    ts[i][j] = -ts[i][j];
    }
  }
}



Resum::Resum(ISR_Handler *const isr,
	     Model_Base *const model):
  Shower_Base("Resum"), p_ampl(NULL)
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

Resum::~Resum()
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

int Resum::PerformShowers()
{
  DEBUG_FUNC(this);
  if (p_ampl==NULL) THROW(fatal_error,"No process info");
  msg_Debugging()<<*p_ampl<<"\n";  
  m_weight=1.0;
  Vec4D_Vector moms(p_ampl->Legs().size());
  Flavour_Vector flavs(p_ampl->Legs().size());
  for (size_t i(0);i<p_ampl->Legs().size();++i) {
    moms[i]=i<p_ampl->NIn()?-p_ampl->Leg(i)->Mom():p_ampl->Leg(i)->Mom();
    flavs[i]=i<p_ampl->NIn()?p_ampl->Leg(i)->Flav().Bar():p_ampl->Leg(i)->Flav();
  }
  m_rn[0]=ran->Get();
  m_rn[1]=ran->Get();
  for (size_t n(0);n<m_obss.size();++n) {
    DEBUG_FUNC(m_obss[n]->Name());
    m_a.clear();
    m_b.clear();
    m_logdbar.clear();
    for (size_t i(0);i<moms.size();++i) {
      Obs_Params ps=m_obss[n]->Parameters
	(&moms.front(),&flavs.front(),moms.size(),i);
      m_a.push_back(ps.m_a);
      m_b.push_back(ps.m_b);
      m_logdbar.push_back(ps.m_logdbar);
    }
    size_t i=1+m_hist[n]->Nbin()*ran->Get();
    double xl=m_hist[n]->LowEdge(i), yl=Value(xl,n);
    double xh=m_hist[n]->HighEdge(i), yh=Value(xh,n);
    m_ress[n][0]=i+1;
    m_ress[n][1]=m_weight*(yh-yl)*m_hist[n]->Nbin();
    msg_Debugging()<<"Bin["<<i<<"]("<<xl<<","<<xh<<"): "<<yl<<" "<<yh<<"\n";
  }  
  CleanUp();
  return 1;
}

double Resum::Value(const double &v, const int n)
{
  DEBUG_FUNC(v);
  double L=log(1.0/v);
  double Rp=0.0, Collexp=0.0, Softexp=0.0, PDFexp=0.0;
  double weight=CalcS(L,Softexp);
  //weight*=1 //non-global logs  
  //weight*=(1+delta) //finite aS corrections
  //calc PDF factor for IS legs
  weight*=CalcPDF(L, PDFexp); 
  //calc collinear piece
  weight*=exp(CalcColl(L,1,Rp,Collexp)); 
  weight*=m_obss[n]->CalcF(Rp);
  if (m_amode) {
    weight=0.0;
    if (m_amode&1) weight+=Collexp+Softexp;
    if (m_amode&2) weight+=PDFexp;
    return weight;
  }
  return weight;
}

int Resum::PerformDecayShowers()
{
  DEBUG_FUNC(this);
  return PerformShowers();
}

bool Resum::ExtractPartons(Blob_List *const bl)
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

void Resum::CleanUp()
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

Cluster_Definitions_Base *Resum::GetClusterDefinitions()
{
  return NULL;
}

bool Resum::PrepareShower
(Cluster_Amplitude* ampl,const bool &soft)
{
  DEBUG_FUNC(this);
  DEBUG_VAR(ampl->Proc<Process_Base>());

  p_ampl=ampl->Copy();
  p_ampl->SetNIn(0);
  
  //Compute invariant s12 for pre-flavor sorted process
  s_12 = sqrt((p_ampl->Leg(0)->Mom()+p_ampl->Leg(1)->Mom())*(p_ampl->Leg(0)->Mom()+p_ampl->Leg(1)->Mom())); 


  std::string pname=Process_Base::GenerateName(p_ampl);
  Process_Base::SortFlavours(p_ampl);
  CMetric_Base* cmetric;

  n_g=0;
  n_q=0;
  n_aq=0;
   
  for (size_t i=0;i<p_ampl->Legs().size();i++) {
    Flavour flav = p_ampl->Leg(i)->Flav();
    if (i<2) flav=flav.Bar();
    if (flav==Flavour(kf_gluon)) n_g++;
    if (flav.IsQuark() && !flav.IsAnti()) n_q++;
    if (flav.IsQuark() && flav.IsAnti()) n_aq++;
  }

  //number of color singlets
  color_sings = p_ampl->Legs().size() - (n_g + n_aq + n_q);
   
  std::string name= ToString(n_g)+"_G_"+ToString(n_q)+"_Q_"+ToString(n_aq)+"_AQ";
  
  msg_Debugging()<<" Found process with "<<n_g<<" Gluons, "<<n_q<<" Quarks and "<<n_aq<<" Anti-Quarks : "<<name<<std::endl;

  //search for metric in m_cmetrics, 
  //otherwise get it from available getters
  CMetric_Map::const_iterator CMiter=m_cmetrics.find(pname);
  //get_Cbasis(name)
  if (CMiter!=m_cmetrics.end()) {
    cmetric=CMiter->second;
    msg_Debugging()<<" found metric in list : "<<pname<<" -> "<<cmetric<<std::endl;
  }
  else {
    //initial bases calc 
    //Compute metric for process arranged like q...qb...g
    msg_Debugging() << name << std::endl;
    cmetric=CMetric_Base::GetCM(CMetric_Key(name,p_ampl));
    if (cmetric==NULL) THROW(not_implemented,"No metric for "+name);
    msg_Debugging()<<"Metric for '"<<pname<<"' is "<<cmetric<<"\n";
    m_cmetrics.insert(make_pair(pname,cmetric));
  }
  
  p_cmetric=cmetric;
 
  //extract colored kinematics
  std::vector<Vec4D> lab_moms;
  std::vector<Flavour> tmp_lab;
  std::vector<int> tmp_sig;
  int p_it;
  for (int i = 0 ; i < (n_g + n_aq + n_q) ; i++) {
    p_it = p_cmetric->Map(i);   
    Flavour flav = p_ampl->Leg(p_it + color_sings)->Flav();
    Vec4D tmp = p_ampl->Leg(p_it + color_sings)->Mom();
    //Determine color sign
    int sign = (tmp[0] > 0 ? 1 : -1);
    tmp_sig.push_back(sign);
    msg_Debugging()<< p_it <<": " << tmp[0] << "  "<<  p_ampl->Leg(p_it)->Flav() << "  " << sign << std::endl;
    if(tmp[0] < 0) tmp=-tmp;
    lab_moms.push_back(tmp);
    tmp_lab.push_back(flav);
  }
  
  Poincare cms_boost(lab_moms[0]+lab_moms[1]);
  
  for (size_t i=0;i<lab_moms.size();i++) {
    Vec4D tmp=lab_moms[i];
    cms_boost.Boost(tmp);
    m_cms.push_back(tmp);
  }

  //kinematic invariants
  for (size_t i=0;i<lab_moms.size();i++) {
      for (size_t j=i+1;j<lab_moms.size();j++) {
	m_Qij.push_back(sqrt((lab_moms[i]+lab_moms[j])*(lab_moms[i]+lab_moms[j])));
	momlabels.push_back(lab_moms[i]);
	momlabels.push_back(lab_moms[j]);
	flavlabels.push_back(tmp_lab[i]);  
	flavlabels.push_back(tmp_lab[j]);
	signlabels.push_back(tmp_sig[i]);  
	signlabels.push_back(tmp_sig[j]);
      }
  }

  p_ampl->Delete();
  p_ampl=ampl->Copy();

  return true;
} 
  
size_t Resum::AddObservable(Observable_Base *const obs,
			    ATOOLS::Histogram *const h)
{
  m_obss.push_back(obs);
  m_hist.push_back(h);
  m_ress.push_back(std::vector<double>(2));
  return m_ress.size()-1;
}

double Resum::CplFac(const ATOOLS::Flavour &fli,const ATOOLS::Flavour &flj,
		    const ATOOLS::Flavour &flk,const int type,
		    const int cpl,const double &mu2) const
{
  THROW(not_implemented,"");
  return -1.0;
}


double Resum::T(const double x)
{
  return -1./p_as->Beta0(p_ampl->MuR2())*log(1.-2.*x);
} 



double Resum::CalcS(const double L, double &Softexp) 
{

  //Exception for n_colored = 2
  if(n_g + n_q + n_aq == 2) return 1.;

  double as = (*p_as)(p_ampl->MuR2());
  double beta0 = p_as->Beta0(p_ampl->MuR2())/M_PI;
  double lambda=as*beta0*L; 
  double t = T(lambda/m_a[0]);
  //order a_s expansion of t
  if (m_amode&1) t = 2*as*L/M_PI;
  size_t numlegs = n_g + n_q + n_aq;

  Hard_Matrix *h(p_comix->ComputeHardMatrix(p_ampl,p_cmetric->Perms()));
  if (h==NULL) msg_Error()<<METHOD<<"(): Error computing H."<<std::endl;
  
  //Color Metric
  msg_Debugging()<<"color metric yields : "<<std::endl;
  std::vector< std::vector< double > > met = p_cmetric->CMetric(), ICmetric = p_cmetric->Imetric() ;
  printMat(met);
  size_t dim = met.size();

  //Inverse metric
  msg_Debugging() << "Inverse Metric" << std::endl;
  printMat(ICmetric);
  
  //Check metric.Inverse = identity
  double value;
  msg_Debugging() << "Check metric.Inverse = identity" << std::endl;
  for(size_t i = 0; i<met.size();i++){
    for(size_t j = 0; j<met.size();j++){
      value = 0.;
      for(size_t k=0; k<met.size();k++){
	value = value + met[i][k]*ICmetric[k][j];
      }
      msg_Debugging() << ( (fabs(value) > .000001 ) ? value : 0 ) << " "; 
    }
    msg_Debugging() << std::endl;
  }
  msg_Debugging() << std::endl;
  
 
  //get the t products & calc Gamma
  std::vector< std::vector< std::vector< double > > > Tprods = p_cmetric->Tprods(); ;
  
  //Check color conservation
  double check;
  double Casmir;
  size_t k_t,i_t, j_t, loop,count;
  for(size_t k = 0; k<numlegs; k++){ 
    msg_Debugging() << "Should be T" << k+1 << ".J = 0" << std::endl;
    k_t = (k < 2 ? k : 2*k-1);
    Casmir = flavlabels[k_t]==Flavour(kf_gluon) ? s_CA : s_CF;
    msg_Debugging() << "is a " << flavlabels[k_t] << std::endl;
    for(size_t i = 0; i<dim; i++){
      for(size_t j = 0; j<dim; j++){
	check = Casmir*met[i][j];
        msg_Debugging()<<"Casimir*met = "<<Casmir<<"*"<<met[i][j]<<" = "<<check<<"\n";
        i_t=0; j_t=0; loop=0; count = 0;
	for(size_t l = 0; l<Tprods.size(); l++) {
	  i_t = loop;
	  j_t = count+loop+1;
	  count++;
	  if(j_t == numlegs-1){ loop++; count = 0; }
	  if(i_t == k || j_t == k){
            msg_Debugging()<<"check + T = "<<check<< " + " << Tprods[l][i][j];
            check += Tprods[l][i][j];
            msg_Debugging()<<" = "<<check<<"\n";
          }

	}
	msg_Debugging() << (fabs(check) > .0001 ? check : 0) << " ";
      }
      msg_Debugging() << std::endl;
    }
    msg_Debugging() << std::endl;
  }
  msg_Debugging() << std::endl;

  //Build Gamma
  Complex *Gamma = new Complex[dim*dim];
  if (Tprods.size()>0) {
    Complex entry;
    msg_Debugging()<<"Gamma = \n";
    for (size_t i=0;i<dim;i++) {
      msg_Debugging()<<"  {";
      for (size_t j=0;j<dim;j++) {
      	  entry = 0;   
	  //Sum over dipoles
          for (size_t k = 0; k < Tprods.size(); k++){
	    entry += - 2.* log(m_Qij[k]/s_12) * Tprods[k][i][j];
	    
	    //Coulomb Gluons 
	    if(signlabels[2*k]*signlabels[2*k+1] == 1) entry += Complex(0.,M_PI)*(Tprods[k][i][j]);
	  }
         
	  Gamma[i*dim+j] = -(t/2.)*entry;
	  msg_Debugging()<<Gamma[i*dim+j].real()<<"+"
	 	         <<Gamma[i*dim+j].imag()<<((j+1<dim)?"I,":"I");
	  }
      msg_Debugging()<<((i+1==dim)?"}\n":"},\n");
    }
    msg_Debugging()<<"}\n";
    }

    //exponentiate gamma with inverse metric
    Complex *GammaI = new Complex[dim*dim];
    Complex *eGamma;
    Complex *cGamma = new Complex[dim*dim];
    matMult(GammaI,Gamma,ICmetric);
    eGamma = c8mat_expm1 (dim,GammaI);
    Complex *hconj=H_conjugate(eGamma,dim);
    matMult(cGamma,hconj,met);
    
    //Soft matrix 
    std::vector< std::vector< double > > Soft(dim);
    for (size_t i=0;i<dim;i++) {
      Soft[i].resize(dim);
      for (size_t j=0;j<dim;j++) {
	for (size_t k=0;k<dim;k++) {
	  Soft[i][j] = Soft[i][j] + real(cGamma[i*dim+k]*eGamma[k*dim+j]);  
	}
      }
    }
    msg_Debugging()<<"Soft  \n";	
    
    //Hard Matrix in T-basis
    std::vector< std::vector< double > > Hard(dim);
    for(size_t i = 0; i<dim; i++){  
        Hard[i].resize(dim);
	for( size_t j = 0; j<dim; j++ ){
	  Hard[i][j] = (i < p_cmetric->Perms().size() && j < p_cmetric->Perms().size()) ? access(*h,i,j) : 0; 
	}
    }
    msg_Debugging()<<"Hard\n";
    printMat(Hard);    

    //Trace of hard matrix (Hard Matrix Element)
    //Note that normalization of H drops out in S
    double traceH = 0;
    for(size_t i = 0; i<dim; i++){
	for(size_t j = 0; j<dim; j++){
	traceH = traceH + Hard[j][i]*met[i][j];
        }
    }
    
    //Leading order expansion Tr(-t*H*Gamma);
    double traceHG = 0;
    for(size_t i = 0; i<dim; i++){
	for(size_t j = 0; j<dim; j++){
	traceHG = traceHG + 2*Hard[j][i]*real(Gamma[i*dim+j]);
	}
    }
  
    Softexp=traceHG/traceH/m_a[0];

    //Hard-soft contraction 
    double traceHS = 0;
    for(size_t i = 0 ; i<dim ; i++ ){
	for( size_t j = 0 ; j<dim ; j++ ){
	traceHS = traceHS + Hard[j][i]*Soft[i][j];
        }
    }

    msg_Debugging()<<"==============================================" << std::endl;
    msg_Debugging()<<"Soft function:" << std::endl;
    msg_Debugging()<<"==============================================" << std::endl;
    
    msg_Debugging()<< "alpha_s: " << as << std::endl;
    msg_Debugging()<< "evolution variable t: " << t << std::endl;
    msg_Debugging()<< "Log(1/v): " << L << std::endl;
    msg_Debugging()<< std::endl;
    msg_Debugging()<< "Kinematics" << std::endl;
    msg_Debugging()<< *p_ampl << std::endl;
    
    msg_Debugging()<< "Check energy-momentum conservation: ";
    Vec4D testEM(0,0,0,0);
    for(size_t i = 0; i < numlegs; i++){
      testEM = testEM + p_ampl->Leg(i)->Mom();
    }
    msg_Debugging()<< testEM << std::endl;
    msg_Debugging()<<"Tr( c H ): " << traceH << std::endl;
    msg_Debugging()<<"Softexp: " << Softexp << std::endl;
    msg_Debugging()<<"Tr( H G Gb ) / Tr( c H ): " << traceHS/traceH << std::endl;
    msg_Debugging() << std::endl;

     //Print Tprods
     size_t coin = 0;
     for(size_t i = 0; i<numlegs;i++){
     for(size_t j = i+1; j<numlegs;j++){
       double traceHT=0;
       //Calculate trace with Hard
       double traceHS = 0;
       for(size_t l = 0 ; l<dim ; l++ ){
   	 for( size_t m = 0 ; m<dim ; m++ ){
	   traceHT = traceHT + Hard[m][l]*Tprods[coin][l][m];
           }
       }

     msg_Debugging() << "T" << i << ".T" << j << "  :  " << std::endl;
     msg_Debugging() << "T" << i << " Flavour: " << flavlabels[2*coin] << ":  four vec: " << momlabels[2*coin] << std::endl;
     msg_Debugging() << "T" << j << " Flavour: " << flavlabels[2*coin+1] << ":  four vec: " << momlabels[2*coin+1] << std::endl;
     msg_Debugging() << "log(Qij/Q12): " << log(m_Qij[coin]/s_12) << std::endl;
     msg_Debugging() << "Qij*Qij: " << std::setprecision(9) << m_Qij[coin]*m_Qij[coin] << std::endl;
     msg_Debugging() <<  "Tr( T.H )/Tr(c.H): " << traceHT/traceH << std::endl;
     printMat(Tprods[coin]);
     msg_Debugging() << std::endl;
     coin++;
      }
    }

    msg_Debugging()<<"==============================================" << std::endl;
    msg_Debugging()<<"end checks" << std::endl;
    msg_Debugging()<<"==============================================" << std::endl;
    msg_Debugging()<< std::endl;
 
    delete [] eGamma;
    delete [] hconj;
    delete [] Gamma;
    delete [] GammaI;
    delete [] cGamma;
    p_comix->Reset();
    
    return traceHS/traceH;
  
}

double Resum::CalcF(const double Rp) 
{
  //implements F function for thrust 
  //use implementation of Ln(Gamma) in MathTools
  double F = exp(-GAMMA_E*Rp-Gammln(1.+Rp));
  return F;
}

double Resum::CalcColl(const double L,const int order,double &Rp, double &Collexp) 
{

  const double nf = double(p_as->Nf(p_ampl->MuR2()));
  const double Tf = s_TR*nf;
  const double beta0 = p_as->Beta0(p_ampl->MuR2())/M_PI;
  const double beta1 = p_as->Beta1(p_ampl->MuR2())/(M_PI*M_PI);
  const double K_CMW = s_CA*(67./18. - M_PI*M_PI/6.)- Tf*10./9.;
  
  double R=0;

  Poincare cms(p_ampl->Leg(0)->Mom()+p_ampl->Leg(1)->Mom());
  
  for(size_t i =0;i<p_ampl->Legs().size();i++) 
    {
    
      double colfac = 0.;
      double hardcoll=0.;
      double as;
      double El;

      as = (*p_as)(p_ampl->MuR2());
      Vec4D pl(p_ampl->Leg(i)->Mom());
      cms.Boost(pl);
      El = dabs(pl[0]);

      if (p_ampl->Leg(i)->Flav().StrongCharge()==8) {colfac = s_CA; hardcoll=-M_PI*beta0/s_CA;}
      else if (abs(p_ampl->Leg(i)->Flav().StrongCharge())==3) {colfac = s_CF; hardcoll=-3./4.;} 
      else continue;

      double t_scale = 0.5*(sqrt(pow(p_ampl->Leg(2)->Mom()[1],2) + pow(p_ampl->Leg(2)->Mom()[2],2)+
				 pow(p_ampl->Leg(3)->Mom()[1],2) + pow(p_ampl->Leg(3)->Mom()[2],2)));
      double Q=sqrt(p_ampl->MuQ2());
      double Q12=s_12;
      double lambda=as*beta0*L;
      double Lmur=log((p_ampl->MuR2())/(Q*Q));

      //The following formulae are taken from Appendix A of hep-ph/0407286. 

      if (abs(m_b[i])>0.0001)
	{
	  if (order>=0) {    
	    //LL part
	    double r1= 1./2./M_PI/pow(beta0,2.)/as/m_b[i]*((m_a[i]-2.*lambda)*log(1.-2.*lambda/m_a[i])
							    -(m_a[i]+m_b[i]-2.*lambda)*log(1.-2.*lambda/(m_a[i]+m_b[i])));
	    R+=((-1.)*colfac*r1);
	  }
	  if (order>=1) {	    
	  //NLL part   //note Lmur term in r2_cmw
	    double r2_cmw=(K_CMW/pow(2.*M_PI*beta0,2.)-Lmur/M_PI/beta0/2.)*((m_a[i]+m_b[i])*log(1.-2.*lambda/(m_a[i]+m_b[i]))
						       -m_a[i]*log(1.-2.*lambda/m_a[i]));
	    double r2_beta1=beta1/2./M_PI/pow(beta0,3.)*(m_a[i]/2.*pow(log(1-2.*lambda/m_a[i]),2.)
							 -0.5*(m_a[i]+m_b[i])*pow(log(1.-2.*lambda/(m_a[i]+m_b[i])),2.)
							 +m_a[i]*log(1-2.*lambda/m_a[i])
							 -(m_a[i]+m_b[i])*log(1.-2.*lambda/(m_a[i]+m_b[i])));
	    double r2=1./m_b[i]*(r2_cmw+r2_beta1);
	    double r1p=1./m_b[i]*(T(lambda/m_a[i])-T(lambda/(m_a[i]+m_b[i])));	    
	    R+=(-1.)*colfac*(r2+r1p*(m_logdbar[i]+m_a[i]*log(Q/Q12)-m_b[i]*log(2.0*El/Q))+hardcoll*T(lambda/(m_a[i]+m_b[i]))+log(Q12/Q)*T(lambda/m_a[i]));
	    Rp+=r1p*colfac;
	    
	  }
	}
      else
	{
	  if (order>=0) {    
	    //LL part
	    double r1= -1./2./M_PI/pow(beta0,2.)/as*(2.*lambda/m_a[i]+log(1.-2.*lambda/m_a[i]));
	    R+=(-1.)*colfac*r1;
	  }
	  if (order>=1) {	    
	  //NLL part
	    double r2_cmw=(K_CMW/pow(2.*M_PI*beta0,2.)-Lmur/M_PI/beta0/2.)*(log(1.-2.*lambda/m_a[i])+2./m_a[i]*lambda/(1.-2./m_a[i]*lambda));
	    double r2_beta1=-beta1/2./M_PI/pow(beta0,3.)*(1./2.*pow(log(1-2.*lambda/m_a[i]),2.)
							  +(log(1-2.*lambda/m_a[i])+2./m_a[i]*lambda)/(1.-2*lambda/m_a[i]));
	    double r2=(r2_cmw+r2_beta1);
	    double r1p=2./(m_a[i]*m_a[i])/(M_PI*beta0)*lambda/(1.-2.*lambda/m_a[i]);
	    R+=(-1.)*colfac*(r2+r1p*(m_logdbar[i]+m_a[i]*log(Q/Q12))+hardcoll*T(lambda/m_a[i])+log(Q12/Q)*T(lambda/m_a[i]));
	    Rp+=r1p*colfac;
	  }
	}

      Collexp+= -2./M_PI*as*(colfac) * ( L/2.0/m_a[i]/(m_a[i]+m_b[i]) + hardcoll/(m_a[i]+m_b[i]) + 1./m_a[i]/(m_a[i]+m_b[i])*(m_logdbar[i]+m_a[i]*log(Q/Q12)-m_b[i]*log(2.0*El/Q)) + log(Q12/Q)/m_a[i])*L;
    }
  return R;
}

double Resum::CalcPDF(const double L, double &PDFexp) 
{
  //strong coupling & PDFs
  double as = (*p_as)(p_ampl->MuR2());

  double old_pdffac=1.,new_pdffac = 1.;
  
  double scale= p_ampl->MuF2();
  msg_Debugging() << "scale before: " << scale << std::endl;

  for (size_t i=0;i<2;i++) {
    if (p_ampl->Leg(i)->Flav().IsLepton()) continue;
    double x = 0;

    if (i==0) x=(-p_ampl->Leg(i)->Mom()).PPlus()/rpa->gen.PBeam(0).PPlus();
    if (i==1) x=(-p_ampl->Leg(i)->Mom().PMinus())/rpa->gen.PBeam(1).PMinus();
    
    //original PDF
    p_pdf[i]->Calculate(x,scale);
    //O(as) expansion of the PDF evolution
    Single_Process *proc(p_ampl->Proc<Single_Process>());
    double z(x+(1.0-x)*m_rn[i]);
    PDFexp+=-2.0/(m_a[i]+m_b[i])*proc->CollinearCounterTerms(i,p_ampl->Leg(i)->Flav().Bar(),-p_ampl->Leg(i)->Mom(),z,0.0,0.0,1)*L;

    old_pdffac*=p_pdf[i]->GetXPDF(p_ampl->Leg(i)->Flav().Bar());
    //new PDF scale
    double scalefac = pow(exp(-L),2./(m_a[i]+m_b[i]));
    if (scale*scalefac<p_pdf[i]->Q2Min()) {
      //freeze PDF at Q2Min 
      p_pdf[i]->Calculate(x,p_pdf[i]->Q2Min());
    }
    else {
      p_pdf[i]->Calculate(x,scale*scalefac);
    }
    msg_Debugging() << "scale after: " << scale << std::endl;
    new_pdffac*=p_pdf[i]->GetXPDF(p_ampl->Leg(i)->Flav().Bar());      
  }

  return new_pdffac/old_pdffac;
}

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



DECLARE_GETTER(Resum,"Resum",Shower_Base,Shower_Key);

Shower_Base *Getter<Shower_Base,Shower_Key,Resum>::
operator()(const Shower_Key &key) const
{
  return new Resum(key.p_isr,key.p_model);
}

void Getter<Shower_Base,Shower_Key,Resum>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"The Colorful Resummation"; 
}
