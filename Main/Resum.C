#include "Main/Resum.H"

#include "PHASIC++/Process/Single_Process.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "MODEL/Main/Model_Base.H"

#include "MODEL/Main/Running_AlphaS.H"

#include "Main/Cluster_Definitions.H"

#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Math/Histogram.H"
#include "Analysis/Observable_Base.H"
#include "Tools/StringTools.H"
#include "Math/Matrix.H"

#include <vector>
#include <complex>
#include <algorithm>

using namespace RESUM;
using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;
using namespace MODEL;

using std::string;

using std::vector;
using std::complex;


double Params::beta0(double scale2) const {
  if(!m_largeNC) return p_as->Beta0(scale2)/M_PI;
  else return 11./12./M_PI;
}
double Params::beta1(double scale2) const {
  if(!m_largeNC) return p_as->Beta1(scale2)/(M_PI*M_PI);
  else return 17./24./M_PI/M_PI;
}
double Params::K_CMW(double scale2) const {
  if(!m_largeNC) return  s_CA*(67./18. - M_PI*M_PI/6.)- m_TR*p_as->Nf(scale2)*10./9.;
  else return 67./18. - M_PI*M_PI/6.;
}
double Params::CollDimGlue(double scale2) const {
  if(!m_largeNC) return -M_PI*beta0(scale2)/s_CA;
  else return -M_PI*beta0(scale2);
}
double Params::CollDimQuark(double scale2) const {
  return -3./4.;
}


Resum::Resum(ISR_Handler *const isr,
	     Model_Base *const model):
  Shower_Base("Resum"), p_ampl(nullptr)
{
  p_clus = new Cluster_Definitions();
  p_as=(Running_AlphaS*)model->GetScalarFunction("alpha_S");
  
  p_pdf = new PDF_Base*[2];
  for (int i=0;i<2; i++) p_pdf[i] = isr->PDF(i);

  Data_Reader read(" ",";","#","=");
  const string& mode = read.GetValue<string>("RESUM::MODE",
                                             read.GetValue<string>("RESUM_MODE","RESUM"));
  if(is_int(mode)) {
    m_amode = static_cast<MODE>(to_type<int>(mode));
  }
  else {
    for(const string& m: split(mode,"\\|")) {
      m_amode = static_cast<MODE>(m_amode | m_ModeToEnum.at(m));
    }
  }
  rpa->gen.SetVariable("SCALES", read.GetValue<string>("SCALES", "VAR{sqr(91.188)}"));
  rpa->gen.SetVariable("RESUM::pre_calc", read.GetValue<string>("RESUM::pre_calc", "pre_calc"));
  
  
  if (rpa->gen.Variable("SHOWER_GENERATOR")=="") {
    rpa->gen.SetVariable("SHOWER_GENERATOR",ToString(this));
  }
  msg_Debugging()<<"Resum Mode: "<<m_amode<<"\n";
  m_params = Params(p_as, (m_amode & MODE::LARGENC));
}

Resum::~Resum()
{
  if (p_ampl) p_ampl->Delete();
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
  if (p_ampl==nullptr) THROW(fatal_error,"No process info");
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
  // loop over observables
  for (size_t n = 0; n<m_obss.size(); n++) {
    DEBUG_FUNC(m_obss[n]->Name());
    m_a.clear();
    m_b.clear();
    m_logdbar.clear();
    for (size_t i(0);i<moms.size();++i) {
      Obs_Params ps = m_obss[n]->Parameters(moms, flavs, i);
      m_a.push_back(ps.m_a);
      m_b.push_back(ps.m_b);

      /// @TODO we have ps.m_d and ps.m_g, in CAESAR notation
      /// logdbar = log (d*average(g))
      /// but we just use logdbar = ps.m_d
      m_logdbar.push_back(ps.m_logdbar);
    }
    m_F = m_obss[n]->FFunction(moms, flavs);
    // select a random bin
    const size_t i = (1+m_hist[n]->Nbin())*ran->Get();    
    const double xl = m_hist[n]->LowEdge(i);
    const double xh = m_hist[n]->HighEdge(i);
    // const double yl = Value(m_obss[n]->LogArg(xl, moms, flavs), -log(xl));
    // const double yh = Value(m_obss[n]->LogArg(xh, moms, flavs), -log(xh));
    m_zcut = m_obss[n]->GroomZcut();
    m_beta = m_obss[n]->GroomBeta();
    m_gmode = m_obss[n]->GroomMode(xl);
    const double ep = m_obss[n]->Endpoint(moms,flavs);
    const double yl = Value(m_obss[n]->LogArg(xl, moms, flavs), -log(m_obss[n]->LogFac()), std::min(xl/ep,1.));
    m_gmode = m_obss[n]->GroomMode(xh);
    const double yh = Value(m_obss[n]->LogArg(xh, moms, flavs), -log(m_obss[n]->LogFac()), std::min(xh/ep,1.));

    // bin to fill
    m_ress[n].first = std::floor(i+1);
    // weight for bin
    m_ress[n].second = m_weight*(yh-yl)*(1+m_hist[n]->Nbin());
    msg_Debugging()<<"Bin["<<i<<"]("<<xl<<","<<xh<<"): "<<yl<<" "<<yh<<" -> "<<(yh-yl)<<"\n";
  }  
  CleanUp();
  return 1;
}


double Resum::Value(const double v, const double LResum, const double epRatio)
{
  DEBUG_FUNC(v);
  if(IsZero(v)) return 0;
  if(v > 1)     return 1;
  const double L = log(1.0/v);
  double Rp = 0.0, CollexpLL=0.0, CollexpNLL=0.0, Softexp=0.0, PDFexp=0.0;
  double weight = 1.;
  weight *= CalcS(L, LResum, Softexp);
  weight *= exp(-epRatio*Softexp);
  // weight*=1 //non-global logs  
  //weight*=(1+delta) //finite aS corrections
  //calc PDF factor for IS legs
  weight *= CalcPDF(L, LResum, PDFexp);
  weight *= exp(-epRatio*PDFexp);
  //calc collinear piece
  weight *= exp(CalcColl(L, LResum, 1, Rp, CollexpLL, CollexpNLL));
  weight *= exp(-epRatio*CollexpNLL);
  if(!std::isnan(Rp)) weight*=m_F(Rp);
  if ((m_amode & (MODE::EXPAND | MODE::PDFEXPAND)) != 0) {
    weight = 0.0;
    if ((m_amode & MODE::COLLEXPAND) != 0) weight += CollexpLL+CollexpNLL*(1.-epRatio);
    if ((m_amode & MODE::SOFTEXPAND) != 0) weight += Softexp*(1.-epRatio);
    if ((m_amode & MODE::PDFEXPAND) != 0)  weight += PDFexp*(1.-epRatio);
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
  if (b==nullptr) THROW(fatal_error,"No Shower blob");
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
  m_comix.Reset();
  if (p_ampl) p_ampl->Delete();
  p_ampl=nullptr;
  
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
  return nullptr;
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
    if (cmetric==nullptr) THROW(not_implemented,"No metric for "+name);
    msg_Debugging()<<"Metric for '"<<pname<<"' is "<<cmetric<<"\n";
    m_cmetrics.insert(make_pair(pname,cmetric));
  }

  /* m_ordered_ids = vector<size_t>(p_ampl->Legs().size()); */
  /* for(size_t i=0; i<p_ampl->Legs().size(); i++) { */
  /*   m_ordered_ids.at(i) = p_ampl->Leg(cmetric->Map(i))->Id(); */
  /* } */
  
  p_cmetric=cmetric;
 
  //extract colored kinematics
  std::vector<Vec4D> lab_moms;
  std::vector<Flavour> tmp_lab;
  std::vector<int> tmp_sig;
  int p_it;
  for (int i = 0 ; i < (n_g + n_aq + n_q) ; i++) {
    p_it = p_cmetric->Map(i);  
    Flavour flav = p_ampl->Leg(p_it)->Flav();
    Vec4D tmp = p_ampl->Leg(p_it)->Mom();
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
  m_ress.push_back({-1,-1});
  return m_ress.size()-1;
}

void Resum::ResetObservables() {
  m_obss.clear();
  m_hist.clear();
  m_ress.clear();
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
  return -1./m_params.beta0(p_ampl->MuR2())/M_PI * log(1.-2.*x);
} 



double Resum::CalcS(const double L, const double LResum, double &Softexp) 
{
  DEBUG_FUNC(L);
  const size_t numlegs = n_g + n_q + n_aq;
  //Exception for n_colored = 2
  if(numlegs == 2) return 1.;

  const double as = (*p_as)(p_ampl->MuR2());
  const double beta0 = m_params.beta0(p_ampl->MuR2());
  const double lambda = as*beta0*L; 
  const double t = T(lambda/m_a[0]);
  const double t_exp = 2*as*L/M_PI;
  
  const MatrixD& met = p_cmetric->CMetric();
  const MatrixC& ICmetric = p_cmetric->Imetric();

  const size_t dim = met.numCols();
  
 
  //get the t products & calc Gamma
  const std::vector<MatrixD>& Tprods = p_cmetric->Tprods();

  if(msg_LevelIsDebugging()) {
    msg_Debugging()<<"Color Metric = "<<std::endl;
    msg_Debugging()<<met<<std::endl;
    msg_Debugging() << "Inverse Metric = " << std::endl;
    msg_Debugging()<<ICmetric<<std::endl;
    
    msg_Debugging() << "Check metric.Inverse = identity" << std::endl;
    msg_Debugging()<<(MatrixC(met)*ICmetric).setFuzzyZeroToZeroInline()<<std::endl;
    
  

    msg_Debugging()<<"Check color conservation: 2 \\sum Tprods + \\sum Casimir*metric = [0].\n";
    MatrixC Tsum(dim, dim);
    for(const MatrixC& T: Tprods) Tsum += T;
    MatrixC Csum(dim, dim);
    for(size_t k=0; k<numlegs; k++) {
      size_t k_t = (k < 2 ? k : 2*k-1);
      double Cl = flavlabels[k_t]==Flavour(kf_gluon) ? m_params.CA() : m_params.CF();
      Csum += Cl*met;
    }
    msg_Debugging()<<(2.*Tsum+Csum).setFuzzyZeroToZeroInline()<<std::endl;
    
    /// @TODO: this was actually a more detailed check
    /* double check; */
    /* double Casmir; */
    /* size_t k_t,i_t, j_t, loop,count; */
    /* for(size_t k = 0; k<numlegs; k++){  */
    /*   msg_Debugging() << "Should be T" << k+1 << ".J = 0" << std::endl; */
    /*   k_t = (k < 2 ? k : 2*k-1); */
    /*   Casmir = flavlabels[k_t]==Flavour(kf_gluon) ? s_CA : s_CF; */
    /*   msg_Debugging() << "is a " << flavlabels[k_t] << std::endl; */
    /*   for(size_t i = 0; i<dim; i++){ */
    /*     for(size_t j = 0; j<dim; j++){ */
    /*       check = Casmir*met[i][j]; */
    /*       msg_Debugging()<<"Casimir*met = "<<Casmir<<"*"<<met[i][j]<<" = "<<check<<"\n"; */
    /*       i_t=0; j_t=0; loop=0; count = 0; */
    /*       for(size_t l = 0; l<Tprods.size(); l++) { */
    /*         i_t = loop; */
    /*         j_t = count+loop+1; */
    /*         count++; */
    /*         if(j_t == numlegs-1){ loop++; count = 0; } */
    /*         if(i_t == k || j_t == k){ */
    /*           msg_Debugging()<<"check + T = "<<check<< " + " << Tprods[l][i][j]; */
    /*           check += Tprods[l][i][j]; */
    /*           msg_Debugging()<<" = "<<check<<"\n"; */
    /*         } */
    
    /*       } */
    /*       msg_Debugging() << (fabs(check) > .0001 ? check : 0) << " "; */
    /*     } */
    /*     msg_Debugging() << std::endl; */
    /*   } */
    /*   msg_Debugging() << std::endl; */
    /* } */
    /* msg_Debugging() << std::endl; */
    
  }
    
  
  //Build Gamma
  MatrixD ReGamma(dim, dim, 0);
  MatrixC Gamma(dim, dim, 0);
  for(size_t k=0; k<Tprods.size(); k++) {
    ReGamma += log(m_Qij[k]/s_12)*Tprods[k];
    if(signlabels[2*k]*signlabels[2*k+1] == 1) {
      Gamma += MatrixC(Tprods[k]);
    }
  }
  MatrixD Gamma_exp = t_exp*ReGamma;
  Gamma *= complex<double>(0,-t/2.*M_PI);
  Gamma += MatrixC(t*ReGamma);
  if(msg_LevelIsDebugging()) {
    msg_Debugging()<<"Gamma = \n"<<Gamma<<"\n";
  }
    
  //Hard Matrix in T-basis
  const MatrixD& Hard = MatrixC(m_comix.ComputeHardMatrix(p_ampl,
                                                   p_cmetric->Perms()),
                         dim, dim, 0).real();
 if(msg_LevelIsDebugging()) {
   msg_Debugging()<<"Hard Matrix =\n"<<Hard.setFuzzyZeroToZero()<<"\n";
 }
  
 //Trace of hard matrix (Hard Matrix Element)
 //Note that normalization of H drops out in S
 const double traceH = Trace(Hard, met);
 //Leading order expansion Tr(-t*H*Gamma);
 const double traceHG = 2.*Trace(Hard, std::move(Gamma_exp));  
 Softexp=traceHG/traceH/m_a[0];

 
 // Calculate Soft matrix
 const MatrixC& eGamma = exp(ICmetric*Gamma.transposeInPlace());
 const MatrixC& cGamma = Conjugate(eGamma)*eGamma;
 const MatrixD& Soft = met*real(std::move(cGamma));
 if(msg_LevelIsDebugging()) {
   msg_Debugging()<<"Soft Matrix = \n"<<Soft.setFuzzyZeroToZero()<<"\n";
 }
 //Hard-soft contraction 
 const double traceHS = Trace(Hard, std::move(Soft));

 if(msg_LevelIsDebugging()) {
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
   for(const auto& l: p_ampl->Legs()) testEM += l->Mom();
   msg_Debugging()<< testEM << std::endl;
   msg_Debugging()<<"Tr( c H ): " << traceH << std::endl;
   msg_Debugging()<<"Softexp: " << Softexp << std::endl;
   msg_Debugging()<<"Tr( H G Gb ) / Tr( c H ): " << traceHS/traceH << std::endl;
   msg_Debugging() << std::endl;
   //Print Tprods
   size_t k = 0;
   for(size_t i = 0; i<numlegs; i++){
     for(size_t j = i+1; j<numlegs; j++){
       const double traceHT = Trace(Tprods[k],Hard);
       msg_Debugging() << "T" << i << ".T" << j << "  :  " << std::endl;
       msg_Debugging() << "T" << i << " Flavour: " << flavlabels[2*k] << ":  four vec: " << momlabels[2*k] << std::endl;
       msg_Debugging() << "T" << j << " Flavour: " << flavlabels[2*k+1] << ":  four vec: " << momlabels[2*k+1] << std::endl;
       msg_Debugging() << "log(Qij/Q12): " << log(m_Qij[k]/s_12) << std::endl;
       msg_Debugging() << "Qij*Qij: " << std::setprecision(9) << m_Qij[k]*m_Qij[k] << std::endl;
       msg_Debugging() <<  "Tr( T.H )/Tr(c.H): " << traceHT/traceH << std::endl;
       msg_Debugging() << std::endl;
       msg_Debugging()<<Tprods[k]<<std::endl;
       msg_Debugging() << std::endl;
       k++;
     }
   }
   msg_Debugging()<<"Going to return traceHS/traceH = "<< traceHS/traceH<<std::endl;
   msg_Debugging() << std::endl;
   msg_Debugging()<<"==============================================" << std::endl;
   msg_Debugging()<<"end checks" << std::endl;
   msg_Debugging()<<"==============================================" << std::endl;
   msg_Debugging()<< std::endl;
 }
 m_comix.Reset();
 return traceHS/traceH;
 
}


double Resum::CalcColl(const double L, const double LResum, const int order, double &Rp, double &CollexpLL, double& CollexpNLL) 
{

  const double muR2 = p_ampl->MuR2();
  const double beta0 = m_params.beta0(muR2);
  const double beta1 = m_params.beta1(muR2);
  const double K_CMW = m_params.K_CMW(muR2);
  
  double R=0;

  Poincare cms(p_ampl->Leg(0)->Mom()+p_ampl->Leg(1)->Mom());
  
  for(size_t i =0; i<p_ampl->Legs().size(); i++) 
    {
    
      double colfac = 0.;
      double hardcoll = 0.;

      const double as = (*p_as)(muR2);
      Vec4D pl(p_ampl->Leg(i)->Mom());
      cms.Boost(pl);
      const double El = dabs(pl[0]);

      if (p_ampl->Leg(i)->Flav().StrongCharge() == 8) {colfac = m_params.CA(); hardcoll=m_params.CollDimGlue(muR2);}
      else if (abs(p_ampl->Leg(i)->Flav().StrongCharge()) == 3) {colfac = m_params.CF(); hardcoll=m_params.CollDimQuark(muR2);} 
      else continue;

      double t_scale = 0.5*(sqrt(pow(p_ampl->Leg(2)->Mom()[1],2) + pow(p_ampl->Leg(2)->Mom()[2],2)+
				 pow(p_ampl->Leg(3)->Mom()[1],2) + pow(p_ampl->Leg(3)->Mom()[2],2)));
      double Q = sqrt(p_ampl->MuQ2());
      double Q12 = s_12;
      double lambda = as*beta0*L;
      double lambdaZ = as*beta0*log(1./m_zcut);
      double lambda2 = as*beta0*log(1./2.);
      double Lmur=log((muR2)/(Q*Q));

      //The following formulae are taken from Appendix A of hep-ph/0407286. 

      if (!IsZero(m_b[i]))
	{
	  if (order>=0) {    
	    //LL part
            double r1 = 0;
            if(m_gmode & GROOM_MODE::SD) {
              r1 = -1./2./M_PI/pow(beta0,2)/as/m_b[i] * (m_b[i]/(1.+m_beta) * (1.-2.*lambdaZ-m_beta*lambda2)*log(1.-2.*lambdaZ-m_beta*lambda2) \
                                                         - (m_a[i]*(1.+m_beta)+m_b[i])/(1.+m_beta) * (1.-2*(1.+m_beta)/(m_a[i]*(1.+m_beta)+m_b[i])*lambda - 2.*m_b[i]/(m_a[i]*(1.+m_beta)+m_b[i])*lambdaZ-2.*m_b[i]/(m_a[i]*(1.+m_beta)+m_b[i])*lambda2)*log(1.-2*(1.+m_beta)/(m_a[i]*(1.+m_beta)+m_b[i])*lambda - 2.*m_b[i]/(m_a[i]*(1.+m_beta)+m_b[i])*lambdaZ-2.*m_b[i]/(m_a[i]*(1.+m_beta)+m_b[i])*lambda2) \
                                                        + (m_a[i]+m_b[i])*(1.-2./(m_a[i]+m_b[i])*lambda)*log(1.-2./(m_a[i]+m_b[i])*lambda) );              
            }
            else {
              r1= 1./2./M_PI/pow(beta0,2.)/as/m_b[i]*((m_a[i]-2.*lambda)*log(1.-2.*lambda/m_a[i])
                                                      -(m_a[i]+m_b[i]-2.*lambda)*log(1.-2.*lambda/(m_a[i]+m_b[i])));
            }
	    R+=((-1.)*colfac*r1);
	  }
	  if (order>=1) {	    
            //NLL part   //note Lmur term in r2_cmw
            double r2_cmw = 0;
            double r2_beta1 = 0;
            double r1p = 0.;
            if(m_gmode & GROOM_MODE::SD) {
              r2_beta1 = -(beta1/4./M_PI/pow(beta0,3.)) * (m_b[i]/(1.+m_beta) * (pow(log(1.-2.*lambdaZ/* -m_beta*lambda2 */),2)+2.*log(1.-2.*lambdaZ/* -m_beta*lambda2 */)) \
                                                           -(m_a[i]*(1.+m_beta)+m_b[i])/(1.+m_beta) * (pow(log(1.-2.*(1.+m_beta)/(m_a[i]*(1.+m_beta)+m_b[i])*lambda - 2.*m_b[i]/(m_a[i]*(1.+m_beta)+m_b[i])*lambdaZ/* -m_beta*m_b[i]/(m_a[i]*m_beta+m_a[i]+m_b[i])*lambda2 */),2)+2.*log(1.-2.*(1.+m_beta)/(m_a[i]*(1.+m_beta)+m_b[i])*lambda - 2.*m_b[i]/(m_a[i]+m_beta+m_a[i]+m_b[i])*lambdaZ-m_beta*m_b[i]/(m_a[i]*m_beta+m_a[i]/* +m_b[i])*lambda2 */)) ) + (m_a[i]+m_b[i])*(pow(log(1.-2./(m_a[i]+m_b[i])*lambda),2)+2.*log(1.-2./(m_a[i]+m_b[i])*lambda)) );
              r2_cmw = K_CMW/pow(2.*M_PI*beta0,2.) * (m_b[i]/(1.+m_beta) * log(1.-2.*lambdaZ/* -m_beta*lambda2 */) \
                                                      - (m_a[i]*(1.+m_beta)+m_b[i])/(1.+m_beta) * log(1.-(2.*(1.+m_beta))/(m_a[i]*(1.+m_beta)+m_b[i])*lambda - 2.*m_b[i]/(m_a[i]*(1.+m_beta)+m_b[i])*lambdaZ /* - m_beta*m_b[i]/(m_a[i]*(1.+m_beta)+m_b[i])*lambda2 */) + (m_a[i]+m_b[i])*log(1.-2./(m_a[i]+m_b[i])*lambda));
              r1p = -(m_a[i]+m_b[i])/(2.*M_PI*m_b[i]*beta0)*(log(1.-(1.+m_beta)*(m_a[i]+m_b[i])/(m_a[i]*(1.+m_beta)+m_b[i])*lambda-2./(m_a[i]*(1.+m_beta)+m_b[i])*lambdaZ/* -m_b[i]*m_beta/(m_a[i]*(1.+m_beta)+m_b[i])*lambda2 */)-log(1-2.*lambda));
            }
            else {
              r2_cmw=(K_CMW/pow(2.*M_PI*beta0,2.)-Lmur/M_PI/beta0/2.)*((m_a[i]+m_b[i])*log(1.-2.*lambda/(m_a[i]+m_b[i]))
                                                                       -m_a[i]*log(1.-2.*lambda/m_a[i]));
              r2_beta1=beta1/2./M_PI/pow(beta0,3.)*(m_a[i]/2.*pow(log(1-2.*lambda/m_a[i]),2.)
                                                    -0.5*(m_a[i]+m_b[i])*pow(log(1.-2.*lambda/(m_a[i]+m_b[i])),2.)
                                                    +m_a[i]*log(1-2.*lambda/m_a[i])
                                                    -(m_a[i]+m_b[i])*log(1.-2.*lambda/(m_a[i]+m_b[i])));
              r1p=1./m_b[i]*(T(lambda/m_a[i])-T(lambda/(m_a[i]+m_b[i])));
            }
            // subtract NLL contribution of scale variation
            double r2_corr = +LResum*r1p;//-(L-LResum)*r1p;
	    double r2=1./m_b[i]*(r2_cmw+r2_beta1+r2_corr);

	    R+=(-1.)*colfac*(r2+r1p*(m_logdbar[i]+m_a[i]*log(Q/Q12)-m_b[i]*log(2.0*El/Q))+hardcoll*T(lambda/(m_a[i]+m_b[i])) + log(Q12/Q)*T(lambda/m_a[i]));
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
	    double r1p=2./(m_a[i]*m_a[i])/(M_PI*beta0)*lambda/(1.-2.*lambda/m_a[i]);
            // subtract NLL contribution of scale variation
            double r2_corr = +LResum*r1p;//-(L-LResum)*r1p;
            // double r2_pow = epRatio * (-2./M_PI*as) * (hardcoll/(m_a[i]+m_b[i]) + 1./m_a[i]/(m_a[i]+m_b[i])*(m_logdbar[i]+m_a[i]*log(Q/Q12)-m_b[i]*log(2.0*El/Q)+LResum) + log(Q12/Q)/m_a[i]) * L;
            double r2=(r2_cmw+r2_beta1+r2_corr);

	    R+= -colfac*(r2+r1p*(m_logdbar[i]+m_a[i]*log(Q/Q12))+hardcoll*T(lambda/m_a[i]) + log(Q12/Q)*T(lambda/m_a[i]));
	    Rp+=r1p*colfac;
	  }
	}
      
      CollexpLL += -2./M_PI*as*(colfac) * L*L/2.0/m_a[i]/(m_a[i]+m_b[i]);
      CollexpNLL += -2./M_PI*as*(colfac) * ((hardcoll/(m_a[i]+m_b[i]) + 1./m_a[i]/(m_a[i]+m_b[i])*(m_logdbar[i]+m_a[i]*log(Q/Q12)-m_b[i]*log(2.0*El/Q)+LResum) + log(Q12/Q)/m_a[i]))*L;

    }
  return R;
}

double Resum::CalcPDF(const double L, const double LResum, double &PDFexp) 
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
