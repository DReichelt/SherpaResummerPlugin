#include "Main/Resum.H"

#include "PHASIC++/Process/Single_Process.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "MODEL/Main/Model_Base.H"

#include "MODEL/Main/Running_AlphaS.H"

#include "Main/Cluster_Definitions.H"

#include "METOOLS/Explicit/NLO_Counter_Terms.H"
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




Resum::Resum(ISR_Handler *const isr,
	     Model_Base *const model):
  Shower_Base("Resum"), p_ampl(nullptr)
{
  p_clus = new Cluster_Definitions();
  p_as=(Running_AlphaS*)model->GetScalarFunction("alpha_S");
  
  p_isr = isr;

  p_pdf = new PDF_Base*[2];
  for (int i=0;i<2; i++) p_pdf[i] = isr->PDF(i);

  Data_Reader read(" ",";","#","=");
  // by default calculate all expansion parts and the resummed distribution
  // also check for RESUM_MODE for backwards compatibility
  const string& mode = read.GetValue<string>("RESUM::MODE",
                                             read.GetValue<string>("RESUM_MODE",
                                                                   "RESUM|EXPAND"));
  if(is_int(mode)) {
    m_amode = static_cast<MODE>(to_type<int>(mode));
  }
  else {
    for(const string& m: split(mode,"\\|")) {
      m_amode = static_cast<MODE>(m_amode | m_ModeToEnum.at(m));
    }
  }
  const string& mmode = read.GetValue<string>("RESUM::MATCHING", "NONE");
  if(is_int(mmode)) {
    m_mmode = static_cast<MATCH_MODE>(to_type<int>(mmode));
  }
  else {
    for(const string& m: split(mmode,"\\|")) {
      m_mmode = static_cast<MATCH_MODE>(m_mmode | m_MModeToEnum.at(m));
    }
  }
  
  rpa->gen.SetVariable("SCALES", read.GetValue<string>("SCALES", "VAR{sqr(91.188)}"));
  rpa->gen.SetVariable("RESUM::pre_calc", read.GetValue<string>("RESUM::pre_calc", "pre_calc"));
  rpa->gen.SetVariable("RESUM::FFUNCTION::VARIATION",read.GetValue<string>("RESUM::FFUNCTION::VARIATION","0"));
  rpa->gen.SetVariable("RESUM::FFUNCTION::PLOT",read.GetValue<string>("RESUM::FFUNCTION::PLOT","0"));
  rpa->gen.SetVariable("RESUM::FFUNCTION::INC",read.GetValue<string>("RESUM::FFUNCTION::INC","1"));
  
  if (rpa->gen.Variable("SHOWER_GENERATOR")=="") {
    rpa->gen.SetVariable("SHOWER_GENERATOR",ToString(this));
  }
  msg_Debugging()<<"Resum Mode: "<<m_amode<<"\n";
  msg_Debugging()<<"Match Mode: "<<m_mmode<<"\n";
  m_params = Params(p_as, (m_amode & MODE::LARGENC));
  m_LogOrd = read.GetValue<size_t>("RESUM::LogOrd", 1);
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
//   DEBUG_FUNC(this);
  if (p_ampl==nullptr) THROW(fatal_error,"No process info");
//   msg_Debugging()<<*p_ampl<<"\n";
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
  for (m_n = 0; m_n<m_obss.size(); m_n++) {
    DEBUG_FUNC(m_obss[m_n]->Name());
    if(!m_init) m_init = true;
    m_a.clear();
    m_b.clear();
    m_logdbar.clear();
    m_etamin.clear();
    for (size_t i=0; i<moms.size(); i++) {
      if(m_obss[m_n] == nullptr)
        THROW(fatal_error,"Observable "+std::to_string(m_n)+" not initalized.");
      Obs_Params ps = m_obss[m_n]->Parameters(p_ampl, i);
      m_a.push_back(ps.m_a);
      m_b.push_back(ps.m_b);

      m_logdbar.push_back(ps.m_logdbar);
      
      m_etamin.push_back(ps.m_etamin);
    }
    
    m_F = m_obss[m_n]->FFunction(moms, flavs, m_params);

    m_zcut = m_obss[m_n]->GroomZcut();
    m_beta = m_obss[m_n]->GroomBeta();

    m_gmode = m_obss[m_n]->GroomMode();
    m_collgmodes = {moms.size(),m_gmode};
    for(size_t i=0; i<m_xvals[m_n].size(); i++) {
      const double x = m_xvals[m_n][i];
      const double ep = m_obss[m_n]->Endpoint(p_ampl);
      const double p = m_obss[m_n]->LogPow(p_ampl);
      if(m_gmode & GROOM_MODE::SD) {
        for(size_t i=0; i<moms.size(); i++) {
          m_collgmodes[i] = m_obss[m_n]->GroomMode(m_obss[m_n]->LogArg(x, p_ampl),
                                                   moms, flavs, i);    
        }
      }
      FillValue(i,m_obss[m_n]->LogArg(x, p_ampl),
                -log(m_obss[m_n]->LogFac(p_ampl)),
                pow(std::min(x/ep,1.),p));
    }
  }
  CleanUp();
  return 1;
}


double Resum::Value(const double v, const double LResum, const double epRatio)
{
  DEBUG_FUNC(v);
  
  
//   Vec4D_Vector moms(p_ampl->Legs().size());
//   Flavour_Vector flavs(p_ampl->Legs().size());
//   for (size_t i(0);i<p_ampl->Legs().size();++i) {
//       moms[i]=i<p_ampl->NIn()?-p_ampl->Leg(i)->Mom():p_ampl->Leg(i)->Mom();
//       flavs[i]=i<p_ampl->NIn()?p_ampl->Leg(i)->Flav().Bar():p_ampl->Leg(i)->Flav();
//   }
  
//   if(IsZero(v)) return 0;
//   if(v > 1)     return 1;
//   const double L = log(1.0/v);
//   double Rp = 0.0, CollexpLL=0.0, CollexpNLL=0.0, Softexp=0.0, PDFexp=0.0;
//   double ExpAtEnd=0.0, RAtEnd=0.0;
//   double Rp0 = 0.0, CollexpLL0=0.0, CollexpNLL0=0.0;
//   double ExpAtEnd0=0.0, RAtEnd0=0.0;
//   double weight = 1.;
  
//   //calc collinear piece
//   weight *= exp(CalcColl(L, LResum, m_LogOrd, Rp, CollexpLL, CollexpNLL, ExpAtEnd, RAtEnd));
//   weight *= exp(-CalcColl(0., LResum, m_LogOrd, Rp0, CollexpLL0, CollexpNLL0, ExpAtEnd0, RAtEnd0));
//   weight *= exp(-epRatio*RAtEnd0*L);
  
//   if(m_LogOrd > 0) {
//     m_gmode = m_obss_n->GroomMode(v, moms, flavs, -1);
// //     std::cout << "1: " << v << " " << CollexpLL << " " << CollexpNLL << " " << Softexp << " " << PDFexp << std::endl;
//     if(!(m_gmode & GROOM_MODE::SD)) {
//       weight *= CalcS(L, LResum, Softexp);
//       weight *= exp(-epRatio*Softexp);
//     }
//     m_gmode = m_obss_n->GroomMode(v, moms, flavs, -2);
//     if(!(m_gmode & GROOM_MODE::SD)) {    
//       //calc PDF factor for IS legs
//       weight *= CalcPDF(L, LResum, PDFexp);
//       weight *= exp(-epRatio*PDFexp);
// //       std::cout << "2: " << v << " " << CollexpLL << " " << CollexpNLL << " " << Softexp << " " << PDFexp << std::endl;
//     }
//     // weight*=1 //non-global logs  
//     //weight*=(1+delta) //finite aS corrections
//     if(!std::isnan(Rp)) weight*=m_F(Rp);
//   }
//   if ((m_amode & (MODE::EXPAND | MODE::PDFEXPAND)) != 0) {
//     weight = 0.0;
//     if ((m_amode & MODE::COLLEXPAND) != 0) weight += CollexpLL+CollexpNLL-epRatio*ExpAtEnd0*L-CollexpLL0-CollexpNLL0;
//     if ((m_amode & MODE::SOFTEXPAND) != 0) weight += Softexp*(1.-epRatio);
//     if ((m_amode & MODE::PDFEXPAND) != 0)  weight += PDFexp*(1.-epRatio);
// =======
//   double Rp = 0.0;
//   MatrixD G(3,3,0);
//   double CollexpLL_LO=0.0;
//   double CollexpNLL_LO=0.0;
//   double CollexpLL_NLO=0.0;
//   double CollexpNLL_NLO=0.0;
//   double SoftexpNLL_LO=0.0;
//   double SoftexpNLL_NLO=0.0;
//   double PDFexp=0.0;
//   double FexpNLL_NLO = 0.0;
//   double weight = 1.;
//   weight *= exp(CalcColl(L, LResum, 1, Rp, G, SoftexpNLL_LO));
//   // double calcs1 = weight*CalcS(L, LResum, SoftexpNLL_LO, SoftexpNLL_NLO,true);
//   weight *= CalcS(L, LResum, SoftexpNLL_LO, SoftexpNLL_NLO);
//   // msg_Out()<<"Check independence of colour inversion: "<<weight-calcs1<<"\n";
//   // if(!IsZero(weight-calcs1,1e-2)) THROW(fatal_error, "Color inversion check failed -> "+ToString(weight-calcs1))
//   // weight*=1 //non-global logs  
//   //weight*=(1+delta) //finite aS corrections
//   //calc PDF factor for IS legs
//   weight *= CalcPDF(L, LResum, PDFexp);
//   //calc collinear piece
//   const double as = (*p_as)(p_ampl->MuR2())/(2.*M_PI);
//   if(m_mmode & MATCH_MODE::ADD) {
//         weight *= exp(-epRatio*SoftexpNLL_LO);
//         weight *= exp(-epRatio*PDFexp);
//         weight *= exp(-epRatio*pow(as,1)*pow(L,1) * ( G(0,0) ));
//   }
//   if(!std::isnan(Rp)) weight*=m_F(Rp,FexpNLL_NLO);
//   else m_F(0,FexpNLL_NLO); // even if Rp diverges its leading order expansion can be determined
//   // msg_Out()<<v<<" "<<weight<<" ";
//   if ((m_amode & (MODE::EXPAND | MODE::PDFEXPAND)) != 0) {
//     weight = 0.0;
//     MatrixD H(4,4,0);
//     if(m_mmode & MATCH_MODE::LO|MATCH_MODE::NLO) {
//       H(0,1) = pow(as,1)*pow(L,2) * ( G(0,1) );
//       H(0,0) = pow(as,1)*pow(L,1) * ( G(0,0) );
//     }
//     if(m_mmode & MATCH_MODE::NLO) {

//       H(1,3) = pow(as,2)*pow(L,4) * ( 0.5*pow(G(0,1),2) );
//       H(1,2) = pow(as,2)*pow(L,3) * ( G(1,2) + G(0,1)*(G(0,0)) );
//       H(1,1) = pow(as,2)*pow(L,2) * ( 0.5*pow(G(0,0),2) + G(1,1) + 4.*FexpNLL_NLO*pow(G(0,1),2)
//                                       );
//     }
//     weight = H.data().sum(); 
//     if(m_mmode & MATCH_MODE::ADD) {
//       weight -= (epRatio*pow(as,1)*pow(L,1) * ( G(0,0) ) + epRatio*SoftexpNLL_LO + epRatio*PDFexp);
//     }
// >>>>>>> NewObservables
  // }
  // return weight;
  THROW(not_implemented, "This function should not be called currently.");
  return 1;
}

void Resum::FillValue(size_t i, const double v, const double LResum, const double epRatio) {
  DEBUG_FUNC(v);
  if(IsZero(v)) return;
  if(v > 1)     return;
  const double L = log(1.0/v);
  double Rp = 0.0;
  // expansion coefficients for collinear part
  MatrixD G(4,4,0);
  MatrixD G0(4,4,0);
  double ExpAtEnd0 = 0.;
  double RAtEnd0 = 0.;
  // expansion coefficents for soft function
  double SoftexpNLL_LO=0.0;
  double SoftexpNLL_NLO=0.0;
  // expansion coefficient for pdf
  double PDFexp=0.0;
  // coefficient F_2
  double FexpNLL_NLO = 0.0;
  // weight for cumulative distribution
  double weight = 1.;

  if(!(m_amode&MODE::IGNCOLL)) {
    msg_Debugging()<<"Calculate collinear piece.\n";
    weight *= exp(CalcColl(L, LResum, 1, Rp, G, SoftexpNLL_LO));
    weight *= exp(-CalcColl(0, LResum, 1, G0, ExpAtEnd0, RAtEnd0));
  }
  msg_Debugging()<<"Weight after coll = "<<weight<<" Rp = "<<Rp<<" L =  "<<L
                 <<" v = "<<m_xvals[m_n][i]<<" bar{d} = "<<exp(m_logdbar[1])<<" "<<exp(m_logdbar[2])<<".\n";
  if(weight > 1e100) {
    msg_Debugging()<<"Large weight in Resum:";
    msg_Debugging()<<*p_ampl<<"\n\n";            
    if(i>0) msg_Debugging()<<"Last Weight = "<<m_resNLL[m_n][i-1]<<"\n";
  }

  if(!(m_amode&MODE::IGNSOFT)) {
    // some checks for colour calculation
    double dummy1, dummy2;
    const double calcs1 = m_amode & MODE::CKINV ?
      weight*CalcS(L, LResum, dummy1, dummy2, MODE::CKINV) : 0;
    const double calcs2 = m_amode & MODE::CKCOUL ?
      weight*CalcS(L, LResum, dummy1, dummy2, MODE::CKCOUL) : 0;
    msg_Debugging()<<"Calculate soft function.\n";
    weight *= CalcS(L, LResum, SoftexpNLL_LO, SoftexpNLL_NLO);
    // evaluate tests if needed
    if(m_amode & MODE::CKINV) {
      msg_Debugging()<<"Check inversion -> "<<weight-calcs1<<".\n";
      if(!std::isnan(weight) and
         !IsZero(weight-calcs1,1e-6)) msg_Error()<<"Color inversion check failed -> "
                                                 <<weight-calcs1<<" rel. dev. = "<<(weight-calcs1)/weight<<". \n\n";
    }
    if(m_amode & MODE::CKCOUL) {
      msg_Debugging()<<"Check for coulomb cancellation -> "<<weight<<" - "<<calcs2<<" = "<<weight-calcs2<<".\n";
      if(!IsZero(weight-calcs2,1e-6)) {
        if(p_ampl->Leg(0)->Flav().StrongCharge()==0 or p_ampl->Leg(1)->Flav().StrongCharge()==0) {
          msg_Error()<<"Coulomb phases did not cancel despite singlet initial state -> "
                     <<(weight-calcs2)<<" ratio = "<<(weight-calcs2)/weight<<".\n";
        }
      }
      if(IsZero(weight-calcs2,1e-6) and p_ampl->Leg(0)->Flav().StrongCharge()!=0 and p_ampl->Leg(1)->Flav().StrongCharge()!=0)
        msg_Out()<<"Complete cancellation of coulomb phases with non-singlet II states ->"
                 <<(weight-calcs2)<<" ratio = "<<(weight-calcs2)/weight<<".\n";
    }
  }
  msg_Debugging()<<"Weight after soft = "<<weight<<"\n";

  if(!(m_amode&MODE::IGNPDF)) {
    msg_Debugging()<<"Calculate PDF.\n";
    weight *= CalcPDF(L, LResum, PDFexp);
  }
  msg_Debugging()<<"Weight after PDF = "<<weight<<".\n";

  // evaluate possible endpoint corrections
  const double as = (*p_as)(p_ampl->MuR2())/(2.*M_PI);
  const double beta0 = m_params.beta0(p_ampl->MuR2());
  if(m_mmode & MATCH_MODE::ADD or m_mmode & MATCH_MODE::DERIV) { 
    weight *= exp(-epRatio*pow(as,1)*pow(L,1)*4./m_a[0]*SoftexpNLL_LO);
    weight *= exp(-epRatio*pow(as,1)*pow(L,1)*PDFexp);
    if(m_mmode & MATCH_MODE::ADD) {
      weight *= exp(-epRatio*pow(as,1)*pow(L,1) * ( G(1,1) ));
    }
    if(m_mmode & MATCH_MODE::DERIV) {
      weight *= exp(-epRatio*pow(L,1)*RAtEnd0);
    }
  }
  msg_Debugging()<<"Weight after ep subtract = "<<weight<<".\n";

  if(!(m_amode&MODE::IGNFFUNC)) {
    msg_Debugging()<<"Calculate FFunction of Rp = "<<Rp<<".\n";
    if(!std::isnan(Rp)) {
      weight*=m_F(Rp,FexpNLL_NLO);
    }
    else m_F(0,FexpNLL_NLO); // if Rp diverges, still get the first expansion coefficient for F
  }
  msg_Debugging()<<"Weight after F = "<<weight<<".\n";

  // store resummed result
  msg_Debugging()<<"Final weight = "<<weight<<"\n";
  m_resNLL[m_n][i] = weight;//std::isnan(weight) ? 0 : weight;

  // calculate expansion
  if(!(m_amode & MODE::COLLEXPAND)) {
    msg_Debugging()<<"Ignore coll. expansion.\n";
    G = MatrixD(4,4,0);
    FexpNLL_NLO = 0;
  }
  if(!(m_amode & MODE::SOFTEXPAND)) {
    msg_Debugging()<<"Ignore soft expansion.\n";
    SoftexpNLL_LO = 0;
    SoftexpNLL_NLO = 0;
  }
  if(!(m_amode & MODE::PDFEXPAND)) {
    msg_Debugging()<<"Ignore pdf expansion.\n";
    PDFexp = 0;
  }

  MatrixD H(4,4,0);
  H(1,2) = pow(as,1)*pow(L,2) * ( G(1,2) );
  H(1,1) = pow(as,1)*pow(L,1) * ( G(1,1) + 4./m_a[0]*SoftexpNLL_LO + PDFexp);
  H(1,0) = pow(as,1)*pow(L,0) * ( G(1,0) );
  double H10 = pow(as,1)*pow(L,0) * ( G0(1,0) );
  m_resExpLO[m_n][i] = H(1,0)+H(1,1)+H(1,2)-H10;
  // TODO: fix
  if(m_mmode & MATCH_MODE::ADD or m_mmode & MATCH_MODE::DERIV) {
    m_resExpLO[m_n][i] -= epRatio*pow(as,1)*pow(L,1) * 4./m_a[0]*SoftexpNLL_LO;
    if(m_mmode & MATCH_MODE::ADD) {
      m_resExpLO[m_n][i] -= epRatio*pow(as,1)*pow(L,1) * G(1,1);
    }
    if(m_mmode & MATCH_MODE::DERIV) {
      m_resExpLO[m_n][i] -= epRatio*pow(as,1)*pow(L,1) * ExpAtEnd0/pow(as,1);
    }
  }  
  // store leading order expansion
  m_resExpLO[m_n][i] = m_resExpLO[m_n][i];//std::isnan(m_resExpLO[m_n][i]) ? 0. : m_resExpLO[m_n][i];
  
  H(2,4) = pow(as,2)*pow(L,4) * ( 0.5*pow(G(1,2),2) );
  H(2,3) = pow(as,2)*pow(L,3) * ( G(2,3) + G(1,2)*(G(1,1) ) );
  // TODO: Fexp handling might be inconsistent with grooming
  H(2,2) = pow(as,2)*pow(L,2) * ( 0.5*pow(G(1,1),2) + G(2,2) +
                                  4.*FexpNLL_NLO*pow(G(1,2),2) );
  m_resExpNLO[m_n][i] = H(2,4) + H(2,3) + H(2,2);
  m_resExpNLO[m_n][i] += pow(as,2)*pow(L,2)*(16./pow(m_a[0],2)*SoftexpNLL_NLO + 4./m_a[0]*SoftexpNLL_LO*(G(1,1)+L*G(1,2)+2.*M_PI*beta0/m_a[0]));

  // TODO: fix
  if(m_mmode & MATCH_MODE::ADD or m_mmode & MATCH_MODE::DERIV) {
    m_resExpNLO[m_n][i] += pow(epRatio*(pow(as,1)*pow(L,1) * (4./m_a[0]*SoftexpNLL_LO+ G(0,0) +PDFexp)),2)/2.;
    if(m_mmode & MATCH_MODE::ADD) {
      m_resExpNLO[m_n][i] += (H(1,1)+H(1,2))*(-epRatio*(pow(as,1)*pow(L,1) * (4./m_a[0]*SoftexpNLL_LO+G(1,1)+PDFexp)));
    }
    if(m_mmode & MATCH_MODE::DERIV) {
      m_resExpNLO[m_n][i] += (H(1,1)+H(1,2))*(-epRatio*(pow(as,1)*pow(L,1) * (4./m_a[0]*SoftexpNLL_LO+ExpAtEnd0/pow(as,1)+PDFexp)));
    }
  }
  // store next-to-leading order expansion
  m_resExpNLO[m_n][i] = m_resExpNLO[m_n][i];//std::isnan(m_resExpNLO[m_n][i]) ? 0. : m_resExpNLO[m_n][i];
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

std::string Resum::AddObservable(const RESUM::Observable_Key& key,
                                 const std::vector<double>& xvals)
{
  DEBUG_FUNC(key.Name());
  Observable_Base* obs = RESUM::Observable_Getter::GetObject(key.Name(),key);
  if(obs != nullptr) {
    m_obss.push_back(obs);
    m_xvals.push_back(xvals);
    m_resNLL.push_back(vector<double>(xvals.size()));
    m_resExpLO.push_back(vector<double>(xvals.size()));
    m_resExpNLO.push_back(vector<double>(xvals.size()));
    m_ress.push_back({-1,-1});
    msg_Debugging()<<"Added "<<obs->Name()<<" as "<<obs->Tag()<<".\n";
  }
  else  msg_Error()<<"Observable not found: "<<key.Name()<<".\n";
  return obs->Tag();
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


const MatrixC& Resum::Gamma() {
  if(m_cacheMatrices) {return m_Gamma;}
  const size_t dim = p_cmetric->CMetric().numCols();

  MatrixD ReGamma(dim, dim, 0);
  m_Gamma = {dim,dim,0};

}

const MatrixD& Resum::Hard() {
  if(m_cacheMatrices) {return m_Hard;}

}



double Resum::CalcS(const double L, const double LResum, double& SoftexpNLL_LO, double& SoftexpNLL_NLO, MODE Check)
{
  DEBUG_FUNC(L);
  if(m_gmode & GROOM_MODE::SD) {
    SoftexpNLL_LO = 0;
    SoftexpNLL_NLO = 0;
    if(m_gmode & GROOM_MODE::SD_SOFT) {
      THROW(not_implemented, "No non-trivial soft function for grooming implemented.");
    }
    msg_Debugging()<<"Ignoring soft function for groomed observables\n";
    return 1.;
  }
  
  const size_t numlegs = n_g + n_q + n_aq;
  //Exception for n_colored = 2
  if(numlegs == 2) return 1.;

  const double as = (*p_as)(p_ampl->MuR2());
  const double beta0 = m_params.beta0(p_ampl->MuR2());
  const double lambda = as*beta0*L; 
  const double t = T(lambda/m_a[0]);
  const double t_exp = 1.;//as*L/M_PI/2.;//2*as*L/M_PI;
  
  const MatrixD& met = p_cmetric->CMetric();
  InverseLowerTriangular(met);
  // TODO: why is this no ref?
  MatrixC ICmetric = p_cmetric->Imetric();
  if(Check & MODE::CKINV) {
    ICmetric += (MatrixC::diagonal(1,ICmetric.numRows())-ICmetric*MatrixC(met))*MatrixC::random(ICmetric.numRows(),ICmetric.numCols());
  }
  
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

    msg_Debugging() << "Check Inverse.metric = identity" << std::endl;
    msg_Debugging()<<(ICmetric*MatrixC(met)).setFuzzyZeroToZeroInline()<<std::endl;

  

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
  MatrixC Gamma_exp = MatrixC(t_exp*ReGamma) \
    - complex<double>(0,-t_exp/2.*M_PI)*Gamma;
  Gamma *= complex<double>(0,-t/2.*M_PI);
  Gamma += MatrixC(t*ReGamma);
  if(msg_LevelIsDebugging()) {
    msg_Debugging()<<"Gamma = \n"<<Gamma<<"\n";
  }
  if(Check & MODE::CKCOUL) {
    Gamma_exp = MatrixC(Gamma_exp.real());
    Gamma = MatrixC(Gamma.real());
  }
    
  //Hard Matrix in T-basis
  MatrixD Hard = MatrixC(m_comix.ComputeHardMatrix(p_ampl,
                                                   p_cmetric->Perms()),
                         dim, dim, 0).real();
  if(p_cmetric->hasTrafo()) {
    const MatrixD& trafo = p_cmetric->TransformationMatrix();
    const MatrixD& trafo_T = Transpose(trafo);
    Hard = trafo*Hard*trafo_T;
  }
  if(msg_LevelIsDebugging()) {
    msg_Debugging()<<"Hard Matrix =\n"<<Hard.setFuzzyZeroToZero()<<"\n";
  }
  
  //Trace of hard matrix (Hard Matrix Element)
  //Note that normalization of H drops out in S
  const double traceH = Trace(Hard, met);
  //Leading order expansion Tr(-t*H*Gamma);
  SoftexpNLL_NLO = SoftexpNLL_LO;
  SoftexpNLL_LO += 2.*Trace(Hard, Gamma_exp.real())/traceH;
  MatrixC conjGamma = Conjugate(Gamma_exp);
  if((m_amode & MODE::SOFTEXPAND) && (m_mmode & MATCH_MODE::NLO)) {
    SoftexpNLL_NLO *= (SoftexpNLL_LO-SoftexpNLL_NLO/2.);
    SoftexpNLL_NLO += 4.*(Trace(Hard,real(conjGamma*ICmetric*conjGamma))
                          + 2.*Trace(Hard,real(conjGamma*ICmetric*Gamma_exp))
                          + Trace(Hard,real(Gamma_exp*ICmetric*Gamma_exp)))/traceH/8.;
  }
 
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
   msg_Debugging()<<"Softexp lo: " << SoftexpNLL_LO << std::endl;
   if((m_amode & MODE::SOFTEXPAND) && (m_mmode & MATCH_MODE::NLO))
     msg_Debugging()<<"Softexp nlo: " << SoftexpNLL_NLO << std::endl;
   msg_Debugging() <<"Tr( H G Gb ) / Tr( c H ): " << traceHS/traceH << std::endl;
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


double Resum::CalcColl(const double L, const double LResum, const int order, double &Rp, MatrixD& G,
                       double& S1, double& ExpAtEnd, double& RAtEnd) 
{
  DEBUG_FUNC(L);
  const double muR2 = p_ampl->MuR2();
  const double beta0 = m_params.beta0(muR2);
  const double beta1 = m_params.beta1(muR2);
  const double K_CMW = m_params.K_CMW(muR2);
  
  double R=0;

  Poincare cms(p_ampl->Leg(0)->Mom()+p_ampl->Leg(1)->Mom());
  
  Vec4D_Vector moms(p_ampl->Legs().size());
  Flavour_Vector flavs(p_ampl->Legs().size());
  for (size_t i(0);i<p_ampl->Legs().size();++i) {
      moms[i]=i<p_ampl->NIn()?-p_ampl->Leg(i)->Mom():p_ampl->Leg(i)->Mom();
      flavs[i]=i<p_ampl->NIn()?p_ampl->Leg(i)->Flav().Bar():p_ampl->Leg(i)->Flav();
  }

  if(m_collgmodes[0] & GROOM_MODE::SD_COLL or
     m_collgmodes[1] & GROOM_MODE::SD_COLL) {
    THROW(not_implemented, "No non-trivial coll. function for groomed initial states implementd");
  }
  for(size_t i = (m_gmode & GROOM_MODE::SD) ? 2 : 0;
      i<p_ampl->Legs().size(); i++) {
    //m_gmode = m_obss_n->GroomMode(exp(-L), moms, flavs, i);
    msg_Debugging()<<"Calculate radiator for leg "<<i<<".\n";
      double colfac = 0.;
      double hardcoll = 0.;

      const double as = (*p_as)(muR2);
      Vec4D pl(p_ampl->Leg(i)->Mom());
      cms.Boost(pl);
      const double El = dabs(pl[0]);
      
      if (p_ampl->Leg(i)->Flav().StrongCharge() == 8) {
        colfac = m_params.CA();
        hardcoll=m_params.CollDimGlue(muR2);
        msg_Debugging()<<"Gluon, Cl = "<<colfac<<" Bl = "<<hardcoll<<".\n";
      }
      else if (abs(p_ampl->Leg(i)->Flav().StrongCharge()) == 3) {
        colfac = m_params.CF();
        hardcoll=m_params.CollDimQuark(muR2);
        msg_Debugging()<<"Quark, Cl = "<<colfac<<" Bl = "<<hardcoll<<".\n";
      } 
      else {
        msg_Debugging()<<"No strong charge, ignoring.\n";
        continue;
      }

      double t_scale = 0.5*(sqrt(pow(p_ampl->Leg(2)->Mom()[1],2) + pow(p_ampl->Leg(2)->Mom()[2],2)+
				 pow(p_ampl->Leg(3)->Mom()[1],2) + pow(p_ampl->Leg(3)->Mom()[2],2)));
      double Q = sqrt(p_ampl->MuQ2());
      double Q12 = s_12;
      
      
      double lambda = as*beta0*L;
      double Lmur=log(muR2/sqr(Q));

      // needed for SD grooming
      double transp = m_obss[m_n]->GroomTransitionPoint(moms, flavs, i);
      double lambdaZ = as*beta0*log(1./transp)/m_a[i];
      double lambda2 = as*beta0*log(1./2.);

      // The following formulae are taken from Appendix A of hep-ph/0407286. 
      if (!IsZero(m_b[i])) {
	  if (order>=0) {    
	    //LL part
            double r1 = 0;
            if(m_gmode & GROOM_MODE::SD and
               m_collgmodes[i] & GROOM_MODE::SD_COLL) {
              r1 = -1./2./M_PI/pow(beta0,2)/as/m_b[i] * (m_b[i]/(1.+m_beta) * (1.-2.*lambdaZ)*log(1.-2.*lambdaZ) \
                                                         - (m_a[i]*(1.+m_beta)+m_b[i])/(1.+m_beta) * (1.-2*(1.+m_beta)/(m_a[i]*(1.+m_beta)+m_b[i])*lambda - 2.*m_b[i]/(m_a[i]*(1.+m_beta)+m_b[i])*lambdaZ)*log(1.-2*(1.+m_beta)/(m_a[i]*(1.+m_beta)+m_b[i])*lambda - 2.*m_b[i]/(m_a[i]*(1.+m_beta)+m_b[i])*lambdaZ) \
                                                         + (m_a[i]+m_b[i])*(1.-2./(m_a[i]+m_b[i])*lambda)*log(1.-2./(m_a[i]+m_b[i])*lambda) );                             
            } // end grooming for LL parts for b != 0
            else {
              msg_Debugging()<<"Argument of log is "<<1.-2.*lambda/m_a[i]<<" and "<<1.-2.*lambda/(m_a[i]+m_b[i])<<".\n";
              r1= 1./2./M_PI/pow(beta0,2.)/as/m_b[i]*((m_a[i]-2.*lambda)*log(1.-2.*lambda/m_a[i])
                                                      -(m_a[i]+m_b[i]-2.*lambda)*log(1.-2.*lambda/(m_a[i]+m_b[i])));
              
            }
            msg_Debugging()<<"LL contribution = "<<-colfac*r1<<"\n";
	    R -= colfac*r1;
	  } // end LL for b != 0
	  if (order>=1) {	    
            //NLL part   //note Lmur term in r2_cmw
            double r2_cmw = 0.;
            double r2_beta1 = 0.;
            double r1p = 0.;
            double r1d = 0.;
            if(m_gmode & GROOM_MODE::SD and
               m_collgmodes[i] & GROOM_MODE::SD_COLL) {
              r2_cmw = (K_CMW/pow(2.*M_PI*beta0,2.)+Lmur/M_PI/beta0/2.) * (m_b[i]/(1.+m_beta) * log(1.-2.*lambdaZ) \
                                                                           - (m_a[i]*(1.+m_beta)+m_b[i])/(1.+m_beta) * log(1.-(2.*(1.+m_beta))/(m_a[i]*(1.+m_beta)+m_b[i])*lambda - 2.*m_b[i]/(m_a[i]*(1.+m_beta)+m_b[i])*lambdaZ) + (m_a[i]+m_b[i])*log(1.-2./(m_a[i]+m_b[i])*lambda));
                
              r2_beta1 = -(beta1/4./M_PI/pow(beta0,3.)) * (m_b[i]/(1.+m_beta) * (pow(log(1.-2.*lambdaZ),2)+2.*log(1.-2.*lambdaZ)) \
                                                           -(m_a[i]*(1.+m_beta)+m_b[i])/(1.+m_beta) * (pow(log(1.-2.*(1.+m_beta)/(m_a[i]*(1.+m_beta)+m_b[i])*lambda - 2.*m_b[i]/(m_a[i]*(1.+m_beta)+m_b[i])*lambdaZ),2)+2.*log(1.-2.*(1.+m_beta)/(m_a[i]*(1.+m_beta)+m_b[i])*lambda - 2.*m_b[i]/(m_a[i]*(1.+m_beta)+m_b[i])*lambdaZ)) + (m_a[i]+m_b[i])*(pow(log(1.-2./(m_a[i]+m_b[i])*lambda),2)+2.*log(1.-2./(m_a[i]+m_b[i])*lambda)) );
                
              r1p = -1./(M_PI*m_b[i]*beta0)*(log(1.-(1.+m_beta)*(m_a[i]+m_b[i])/(m_a[i]*(1.+m_beta)+m_b[i])*lambda-2.*m_b[i]/(m_a[i]*(1.+m_beta)+m_b[i])*lambdaZ)-log(1.-2./(m_a[i]+m_b[i])*lambda));
                
              r1d = -1./(M_PI*m_a[i]*beta0)/(m_beta+1.)*(log(1.-(1.+m_beta)*(m_a[i]+m_b[i])/(m_a[i]*(1.+m_beta)+m_b[i])*lambda-2.*m_b[i]/(m_a[i]*(1.+m_beta)+m_b[i])*lambdaZ)-log(1.-2./m_a[i]*lambdaZ));                          
            } // end grooming for NLL parts for b != 0
            else {
              r2_cmw=(K_CMW/pow(2.*M_PI*beta0,2.)+Lmur/M_PI/beta0/2.)*((m_a[i]+m_b[i])*log(1.-2.*lambda/(m_a[i]+m_b[i]))
						       -m_a[i]*log(1.-2.*lambda/m_a[i]));
              r2_beta1 = beta1/2./M_PI/pow(beta0,3.)*(m_a[i]/2.*pow(log(1-2.*lambda/m_a[i]),2.)
                                                      -0.5*(m_a[i]+m_b[i])*pow(log(1.-2.*lambda/(m_a[i]+m_b[i])),2.)
                                                      +m_a[i]*log(1-2.*lambda/m_a[i])
                                                      -(m_a[i]+m_b[i])*log(1.-2.*lambda/(m_a[i]+m_b[i])));
              r1p= 1./m_b[i]*(T(lambda/m_a[i])-T(lambda/(m_a[i]+m_b[i])));
            } // end of NLL parts for b != 0

            // subtract NLL contribution of scale variation
            const double r2_corr = +LResum*r1p;//-(L-LResum)*r1p;
            // TODO: presumably only one version correct? 
            const double r2=1./m_b[i]*(r2_cmw+r2_beta1)+r2_corr; // was in SD
            //double r2=1./m_b[i]*(r2_cmw+r2_beta1+r2_corr); // was in non-SD. This was in the original version of the code, which I believe is a bug

            // add NLL parts to R and Rp
            R -= colfac*(r2+r1p*(m_logdbar[i]-m_b[i]*log(2.0*El/Q))+hardcoll*T(lambda/(m_a[i]+m_b[i])));
            if(m_collgmodes[i] & GROOM_MODE::SD_COLL) {
              R -= colfac*(r1d*m_logdbar[i] + log(Q12/Q)*T(lambdaZ)); // TODO: not sure I understand the last factor. As this is a contribution of soft wide-angle origin this shift lambda/a->lambdaZ
            }
            else {
              R -= colfac*log(Q12/Q)*T(lambda/m_a[i]);
              R += colfac*m_etamin[i]*T(lambda/m_a[i]);
            } 
            Rp+=r1p*colfac;
	  } // end of NLL for b != 0
      } // end of b != 0
      else { // start b == 0           
        if(m_collgmodes[i] & GROOM_MODE::SD_COLL) {
          THROW(not_implemented,"Groomed b = 0 not implimented.");
        }
        else {
          if (order>=0) {    
            //LL part
            msg_Debugging()<<"Argument of log is "<<1.-2.*lambda/m_a[i]<<".\n";
            double r1= -1./2./M_PI/pow(beta0,2.)/as*(2.*lambda/m_a[i]+log(1.-2.*lambda/m_a[i]));
            msg_Debugging()<<"LL contribution = "<<-colfac*r1<<".\n";
            R -= colfac*r1;
          }
          if (order>=1) {	    
            //NLL part
            double r2_cmw=(K_CMW/pow(2.*M_PI*beta0,2.)+Lmur/M_PI/beta0/2.)*(log(1.-2.*lambda/m_a[i])+2./m_a[i]*lambda/(1.-2./m_a[i]*lambda));
            double r2_beta1=-beta1/2./M_PI/pow(beta0,3.)*(1./2.*pow(log(1-2.*lambda/m_a[i]),2.)
                                                          +(log(1-2.*lambda/m_a[i])+2./m_a[i]*lambda)/(1.-2*lambda/m_a[i]));
            double r1p=2./(m_a[i]*m_a[i])/(M_PI*beta0)*lambda/(1.-2.*lambda/m_a[i]);
            // subtract NLL contribution of scale variation
            double r2_corr = +LResum*r1p;
            double r2=(r2_cmw+r2_beta1+r2_corr);
            
            R += -colfac*(r2+r1p*(m_logdbar[i]-m_b[i]*log(2.0*El/Q))+hardcoll*T(lambda/m_a[i]) + log(Q12/Q)*T(lambda/m_a[i]));
            R += colfac*m_etamin[i]*T(lambda/m_a[i]);
            Rp+=r1p*colfac;
          }
        }
      }
      if(m_collgmodes[i] & GROOM_MODE::SD_COLL) {
        // TODO: sd expansion
//         G(1,2) += -2.*colfac*m_a[i]*m_beta;
//         G(1,1) += -4.*colfac*((m_a[i]+m_b[i])*log(1./transp) + hardcoll/(m_a[i]+m_b[i]) + m_beta/(m_a[i]*(1.+m_beta)+m_b[i])/(m_a[i]+m_b[i])*(m_logdbar[i]+m_a[i]*log(Q/Q12)-m_b[i]*log(2.0*El/Q)+LResum)+m_logdbar[i]/m_a[i]/(m_a[i]*(1.+m_beta)+m_b[i]) );
//         G(1,0) += 4.*colfac*( (m_a[i]+m_b[i])*sqr(log(1./transp))/2.0 - (1./m_a[i]/(m_a[i]*(1.+m_beta)+m_b[i])*(m_logdbar[i]+m_a[i]*log(Q/Q12)-m_b[i]*log(2.0*El/Q)+LResum)-m_logdbar[i]/m_a[i]/(m_a[i]*(1.+m_beta)+m_b[i])+ log(Q12/Q)/m_a[i])*log(1./transp) );
        
        // Fixed type and removed contribution removed from R
        G(1,2) += -2.*colfac*m_beta/(m_a[i]+m_b[i])/(m_a[i]*(1.+m_beta)+m_b[i]);
        G(1,1) += -4.*colfac*(log(1./transp)/m_a[i]/(m_a[i]*(1.+m_beta)+m_b[i]) + hardcoll/(m_a[i]+m_b[i]) + m_beta/(m_a[i]*(1.+m_beta)+m_b[i])/(m_a[i]+m_b[i])*(m_logdbar[i]-m_b[i]*log(2.0*El/Q)+LResum)+m_logdbar[i]/m_a[i]/(m_a[i]*(1.+m_beta)+m_b[i]) );
        G(1,0) += 4.*colfac*(sqr(log(1./transp))/2.0/m_a[i]/(m_a[i]*(1.+m_beta)+m_b[i]) - (1./m_a[i]/(m_a[i]*(1.+m_beta)+m_b[i])*(m_logdbar[i]-m_b[i]*log(2.0*El/Q)+LResum)-m_logdbar[i]/m_a[i]/(m_a[i]*(1.+m_beta)+m_b[i])+ log(Q12/Q)/m_a[i])*log(1./transp) );

        ExpAtEnd += -2./M_PI*as*(colfac)*(m_a[i]+m_b[i])*log(1./transp)/m_a[i]/(m_a[i]+m_b[i])/(m_a[i]*(1.+m_beta)+m_b[i]);
        
        RAtEnd += colfac*log(1.-2.*m_b[i]*lambdaZ/(m_a[i]*(1.+m_beta)+m_b[i]))/M_PI/m_b[i]/beta0;

        double r2_beta1AtEnd = -as*beta1/M_PI/beta0/beta0*m_a[i]*(log(1.-2.*m_b[i]*lambdaZ/(m_a[i]*(1.+m_beta)+m_b[i]))+2.*m_b[i]*lambdaZ/(m_a[i]*(1.+m_beta)+m_b[i]))/(1.-2.*m_b[i]*lambdaZ/(m_a[i]*(1.+m_beta)+m_b[i]));
        double r2_cmwAtEnd = 2.*as*beta0*(K_CMW/pow(2.*M_PI*beta0,2.)+Lmur/M_PI/beta0/2.)*(1./(1.-2.*m_b[i]*lambdaZ/(m_a[i]*(1.+m_beta)+m_b[i]))-1.);
        double r2_hardcollAtEnd = 2.*as/M_PI/(m_a[i]+m_b[i]);
        double r1pAtEnd = -2.*as*beta0/m_b[i]*(1./(m_a[i]+m_b[i])-(1.+m_beta)/(m_a[i]*(1.+m_beta)+m_b[i])/(1.-2.*m_b[i]*lambdaZ/(m_a[i]*(1.+m_beta)+m_b[i])));
        double r1dAtEnd = 2.*as/M_PI/m_a[i]/(m_a[i]*(1.+m_beta)+m_b[i])/(1.-2.*m_b[i]*lambdaZ/(m_a[i]*(1.+m_beta)+m_b[i]));
        
        double r2AtEnd=1./m_b[i]*(r2_cmwAtEnd+r2_beta1AtEnd)+LResum*r1pAtEnd;
        
        ExpAtEnd += -2./M_PI*as*(colfac) * (hardcoll/(m_a[i]+m_b[i]) + m_beta/(m_a[i]*(1.+m_beta)+m_b[i])/(m_a[i]+m_b[i])*(m_logdbar[i]+m_a[i]*log(Q/Q12)-m_b[i]*log(2.0*El/Q)+LResum)+m_logdbar[i]/m_a[i]/(m_a[i]*(1.+m_beta)+m_b[i]) );
        
        RAtEnd += (-1.)*colfac*(r2AtEnd+r1pAtEnd*(m_logdbar[i])+r1dAtEnd*m_logdbar[i]+hardcoll*r2_hardcollAtEnd);        
      } // end expansion for grooming
      else {
        
        ExpAtEnd += -2./M_PI*as*(colfac) * (hardcoll/(m_a[i]+m_b[i]) + 1./m_a[i]/(m_a[i]+m_b[i])*(m_logdbar[i]-m_b[i]*log(2.0*El/Q)+LResum) + log(Q12/Q)/m_a[i]);
        
        RAtEnd += -2./M_PI*as*(colfac) * (hardcoll/(m_a[i]+m_b[i]) + 1./m_a[i]/(m_a[i]+m_b[i])*(m_logdbar[i]-m_b[i]*log(2.0*El/Q)+LResum) + log(Q12/Q)/m_a[i]);
        
        
        G(1,2) += -2./m_a[i] * colfac/(m_a[i]+m_b[i]);
        G(1,1) += -colfac*(4.*hardcoll/(m_a[i]+m_b[i]) + 4./(m_a[i]*(m_a[i]+m_b[i]))*(m_logdbar[i]-m_b[i]*log(2.*El/Q)+LResum));
        G(2,3) += -8.*M_PI*beta0/3./pow(m_a[i],2) * colfac * (2.*m_a[i]+m_b[i])/pow(m_a[i]+m_b[i],2);
        G(2,2) += -colfac*(8.*M_PI*beta0 * (hardcoll/pow(m_a[i]+m_b[i],2) + (2.*m_a[i]+m_b[i])/pow(m_a[i]*(m_a[i]+m_b[i]),2)*(m_logdbar[i]-m_b[i]*log(2.*El/Q)+LResum)) \
                           +2.*(K_CMW+M_PI*beta0*2.*Lmur)/m_a[i]/(m_a[i]+m_b[i]));
        G(1,1) += 4.*colfac*m_etamin[i]/m_a[i];
        G(2,2) += 8.*M_PI*beta0*colfac*m_etamin[i]/m_a[i]/m_a[i];
        
        S1 += -colfac*log(Q12/Q);
      } // end expansion without grooming
    } // end loop over legs
  msg_Debugging()<<"Sum of radiators = "<<R<<".\n";
  return R;
}

double Resum::CollinearCounterTerms(const int i, 
                                    const ATOOLS::Flavour &fl,
                                    const ATOOLS::Vec4D &p,
                                    const double &z,
                                    const double muF2) const {
  // determine ct
  double ct = 0.0; 
  double x = p_isr->CalcX(p);
  Flavour jet(kf_jet);
  double fb=p_isr->PDFWeight((1<<(i+1))|8,p,p,muF2,muF2,fl,fl,0);
  if (IsZero(fb)) {
    msg_Tracking()<<METHOD<<"(): Zero xPDF ( f_{"<<fl<<"}("
		  <<x<<","<<sqrt(muF2)<<") = "<<fb<<" ). Skip.\n";
    return 0.0;
  }

  // skip PDF ratio if high-x sanity condition not fullfilled
  if (dabs(fb)<1.0e-4*log(1.0 - x)/log(1.0 - 1.0e-2)){
    msg_Debugging() << "Invalid pdf ratio, ct set to zero." << std::endl;
    return 0.0;
  }

  msg_Debugging()<<"Beam "<<i<<": z = "<<z<<", f_{"<<fl
		 <<"}("<<x<<","<<sqrt(muF2)<<") = "<<fb<<" {\n";
  for (size_t j(0);j<jet.Size();++j) {
    const double Pf = METOOLS::FPab(jet[j],fl,z);
    const double Ps = METOOLS::SPab(jet[j],fl,z);
    if (Pf+Ps==0.0) continue;
    const double Pi = METOOLS::IPab(jet[j],fl,x);
    const double H = METOOLS::Hab(jet[j],fl);
    const double fa=p_isr->PDFWeight
      (1<<(i+1),p/z,p/z,muF2,muF2,jet[j],jet[j],0);
    const double fc=p_isr->PDFWeight
      (1<<(i+1),p,p,muF2,muF2,jet[j],jet[j],0);
    msg_Debugging()<<"  P_{"<<jet[j]<<","<<fl
		   <<"}("<<z<<") = {F="<<Pf<<",S="<<Ps
		   <<",I="<<Pi<<"}, f_{"<<jet[j]<<"}("
		   <<x/z<<","<<sqrt(muF2)<<") = "<<fa
		   <<", f_{"<<jet[j]<<"}("<<x<<","
		   <<sqrt(muF2)<<") = "<<fc<<"\n";
    if (IsZero(fa)||IsZero(fc)) {
      msg_Tracking()<<METHOD<<"(): Zero xPDF. No contrib from "<<j
                    <<". Skip .\n";
    }
    ct += ((fa/z*Pf+(fa/z-fc)*Ps)*(1.0-x)+fc*(H-Pi))/fb;
  }
  msg_Debugging()<<"} -> "<<ct<<"\n";
  return ct;
}



double Resum::CalcPDF(const double L, const double LResum, double &PDFexp) {
  DEBUG_FUNC(L);
  if(m_gmode & GROOM_MODE::SD) {
    PDFexp = 0;
    if(m_gmode & GROOM_MODE::SD_PDF or
       m_collgmodes[0] & GROOM_MODE::SD_PDF or
       m_collgmodes[1] & GROOM_MODE::SD_PDF) {
      THROW(not_implemented, "No non-trivial pdf contribution for grooming implemented.");
    }
    msg_Debugging()<<"Ignoring pdf function for groomed observables.\n";
    return 1.;
  }
  msg_Debugging()<<"Calculate pdf contribution, no grooming assumed.\n";
  //strong coupling & PDFs
  const double as = (*p_as)(p_ampl->MuR2());

  double old_pdffac = 1.;
  double new_pdffac = 1.;

  const double scale= p_ampl->MuF2();
  msg_Debugging() << "scale before: " << scale << "\n";

  for (size_t i=0;i<2;i++) {
    if (p_ampl->Leg(i)->Flav().IsLepton()) continue;

    const double x = i==0 ? (-p_ampl->Leg(i)->Mom()).PPlus()/rpa->gen.PBeam(0).PPlus() : (-p_ampl->Leg(i)->Mom().PMinus())/rpa->gen.PBeam(1).PMinus();

    //original PDF
    p_pdf[i]->Calculate(x,scale);

    //O(as) expansion of the PDF evolution
    // Single_Process *proc(p_ampl->Proc<Single_Process>());
    // if(proc == nullptr) {
    //    THROW(fatal_error,"Internal error in PDF evaluation.\n");
    // }
    const double z = x+(1.0-x)*m_rn[i];

    msg_Debugging()<<"Calculate PDF expansion with z = "<<z<<".\n";
    // PDFexp+=-2.0/(m_a[i]+m_b[i])*proc->CollinearCounterTerms(i,p_ampl->Leg(i)->Flav().Bar(),-p_ampl->Leg(i)->Mom(),z,exp(1.),1.,1.,1.) * (2.*M_PI)/as;
    PDFexp += -2.0/(m_a[i]+m_b[i])*CollinearCounterTerms(i,p_ampl->Leg(i)->Flav().Bar(),-p_ampl->Leg(i)->Mom(),z,scale);
    // msg_Out()<<CollinearCounterTerms(i,p_ampl->Leg(i)->Flav().Bar(),-p_ampl->Leg(i)->Mom(),z,scale)<<" "
    //          <<METOOLS::CollinearCounterTerms(p_ampl->Leg(i)->Flav().Bar(), x, z,as,exp(1.),1.,scale,p_pdf[i])*as/(2.*M_PI)<<"\n";
    // msg_Out()<<x<<" "<<p_isr->CalcX(-p_ampl->Leg(i)->Mom())<<"\n\n";

    const double fb = p_pdf[i]->GetXPDF(p_ampl->Leg(i)->Flav().Bar());
    old_pdffac *= fb;

    // if (dabs(fb/x)<1.0e-4*log(1.0 - x)/log(1.0 - 1.0e-2)){
    //   msg_Debugging() << "Invalid pdf ratio, ct set to zero." << std::endl;
    //   PDFexp += 0;
    // }
    // else {
    //   PDFexp += -2.0/(m_a[i]+m_b[i])*METOOLS::CollinearCounterTerms(p_ampl->Leg(i)->Flav().Bar(), x, z,2*M_PI,exp(1.),1.,scale,p_pdf[i]);
    // }

    //new PDF scale
    const double scalefac = pow(exp(-L),2./(m_a[i]+m_b[i]));
    if (scale*scalefac<p_pdf[i]->Q2Min()) {
      //freeze PDF at Q2Min 
      p_pdf[i]->Calculate(x,p_pdf[i]->Q2Min());
    }
    else {
      p_pdf[i]->Calculate(x,scale*scalefac);
    }
    msg_Debugging() << "scale after: " << scale*scalefac << std::endl;
    new_pdffac *= p_pdf[i]->GetXPDF(p_ampl->Leg(i)->Flav().Bar());      
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
