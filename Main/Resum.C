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
#include "ATOOLS/Math/MathTools.H"
#include "Analysis/Observable_Base.H"
#include "Tools/StringTools.H"
#include "Math/Matrix.H"
#include "Math/HypGeo.H"
#include "Math/DiGamma.H"

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
using std::valarray;
using std::complex;


Resum::Resum(ISR_Handler *const isr,
	     Model_Base *const model):
  Shower_Base("Resum"), p_ampl(nullptr)
{
  p_clus = new Cluster_Definitions();
  
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
  p_as=(Running_AlphaS*)model->GetScalarFunction("alpha_S");
 
  m_params = Params(p_as, isr->PDF(0), isr->PDF(1), (m_amode & MODE::LARGENC),read.GetValue<double>("RESUM::CONSTAS",-1), 1., 1.);
  m_varParams["CENTRAL"] = m_params;
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


void Resum::SetVariationParameters(const std::vector<SHERPA::Variation_Parameters *> * pvec) {
  if(not pvec) return;
  m_varNames.clear(); 
  for(SHERPA::Variation_Parameters* p: *pvec) {
    m_varNames.push_back(p->m_name);
    m_varParams[p->m_name] = Params(p->p_alphas, p->p_pdf1, p->p_pdf2, m_params.largeNC(), m_params.constAlpha(), 
                                    p->m_muR2fac, p->m_muF2fac);
  }
}

void Resum::SetVariationParameters(const SHERPA::Variations* vars) {
  if(not vars) return;
  SetVariationParameters(vars->GetParametersVector());
}



int Resum::PerformShowers()
{
  if(m_obss.size() == 0) return 1;
  if (p_ampl==nullptr) THROW(fatal_error,"No process info");


  if(m_globaltreat & (TREAT::lo|TREAT::nlo)) {
    std::map<std::string, double> varweights = VarWeightsMap();
    for(size_t i=0; i<m_obss.size(); i++) {
      const std::string& tag = m_obss[i]->Tag();
      msg_Debugging()<<"Fixed order for "<<m_obss[i]->Name()<<" "<<tag<<"\n";
      if(m_resultsLO[tag] and m_treat[i]&TREAT::lo) { 
        m_resultsLO[tag]->FillMCB(m_resFO[i], m_weight, &varweights, 1);
        m_resultsLO[tag]->FillMCB(m_curChannels);
      }
      if(m_resultsNLO[tag] and m_treat[i]&TREAT::nlo) {
        m_resultsNLO[tag]->FillMCB(m_resFO[i], m_weight, &varweights, 1);  
        m_resultsNLO[tag]->FillMCB(m_curChannels);
      }
    }
  }

  if(not (m_globaltreat & TREAT::needs_resum)) return 1;

  Vec4D_Vector moms(p_ampl->Legs().size());
  Flavour_Vector flavs(p_ampl->Legs().size());
  for (size_t i(0);i<p_ampl->Legs().size();++i) {
    moms[i]=i<p_ampl->NIn()?-p_ampl->Leg(i)->Mom():p_ampl->Leg(i)->Mom();
    flavs[i]=i<p_ampl->NIn()?p_ampl->Leg(i)->Flav().Bar():p_ampl->Leg(i)->Flav();
  }
  Poincare cms(p_ampl->Leg(0)->Mom()+p_ampl->Leg(1)->Mom());
  m_rn[0]=ran->Get();
  m_rn[1]=ran->Get();

  m_energy = std::vector<double>(nLegs());
  for (size_t i=0; i<nLegs(); i++) {
    Vec4D pl(p_ampl->Leg(i)->Mom());
    cms.Boost(pl);
    m_energy[i] = std::abs(pl[0]);
  }


  // loop over observables
  for (m_n = 0; m_n<m_obss.size(); m_n++) {
    const std::string& tag = m_obss[m_n]->Tag();
    msg_Debugging()<<"Resummation for "<<m_obss[m_n]->Name()<<" "<<tag<<"\n";
    if(not m_resultsNLL[tag] or not (m_treat[m_n]&TREAT::needs_resum)) continue;
    if(m_obss[m_n] == nullptr)
      THROW(fatal_error,"Observable "+std::to_string(m_n)+" not initalized.");
    if(!m_init) m_init = true;
    m_obss[m_n]->Parameters(p_ampl,m_a,m_b,m_logdbar,m_etamin,m_deltad);


    msg_Debugging()<<"Setting F function.\n";
    m_F = m_obss[m_n]->FFunction(moms, flavs, m_params);

    msg_Debugging()<<"Setting Sngl function.\n";
    m_Sngl = m_obss[m_n]->SnglFunction(moms, flavs, m_params);

    m_nx = m_xvals[m_n].size();

    msg_Debugging()<<"Setting grooming parameters.\n";
    m_zcut = m_obss[m_n]->GroomZcut();
    m_beta = m_obss[m_n]->GroomBeta();

    
    m_gmode = m_obss[m_n]->GroomMode();
    
    m_collgmodes_end = {moms.size(),m_gmode};
    m_softgmode_end = m_gmode;
    if(m_gmode & GROOM_MODE::SD) {
        m_softgmode_end = GROOM_MODE::NONE;
        for(size_t j=0; j<moms.size(); j++) {
            m_collgmodes_end[j] = m_obss[m_n]->GroomMode(1., p_ampl, j);
            //        Transition for soft function if any coll-mode is groomed
            if(m_collgmodes_end[j] & GROOM_MODE::SD_COLL) m_softgmode_end = GROOM_MODE::SD_SOFT;
        }
    }
    
    m_allCollgmodes = std::vector<std::valarray<GROOM_MODE>>(moms.size(),std::valarray<GROOM_MODE>(m_gmode,m_nx));
    m_collGroomed = std::vector<std::valarray<bool>>(moms.size(),std::valarray<bool>(false,m_nx));
    m_softGroomed = std::valarray<bool>(false,m_nx);
    m_allSoftgmodes = std::valarray<GROOM_MODE>(m_gmode,m_nx);

    msg_Debugging()<<"Log rescaling etc.\n";
    m_ep = m_obss[m_n]->Endpoint(p_ampl);
    m_p = m_obss[m_n]->LogPow(p_ampl);
    m_logFac = -log(m_obss[m_n]->LogFac(p_ampl));
    m_epRatio = pow(m_xvals[m_n]/m_ep,m_p);

    m_logArg = m_obss[m_n]->LogArg(m_xvals[m_n], p_ampl);

    for(size_t i=0; i<m_nx; i++) {
      const double x = m_xvals[m_n][i];
      if(m_gmode & GROOM_MODE::SD) {
        m_allSoftgmodes[i] = GROOM_MODE::NONE;
        for(size_t j=0; j<moms.size(); j++) {
          m_allCollgmodes[j][i] = m_obss[m_n]->GroomMode(m_logArg[i], p_ampl, j);
          m_collGroomed[j][i] = m_obss[m_n]->GroomMode(m_logArg[i], p_ampl, j) & GROOM_MODE::SD_COLL;
          //        Transition for soft function if any coll-mode is groomed
          if(m_allCollgmodes[j][i] & GROOM_MODE::SD_COLL){
            if(!(m_softgmode_end & GROOM_MODE::SD_SOFT)) {
              m_softGroomed[i] = true;
              m_allSoftgmodes[i] = GROOM_MODE::SD_SOFT;
            }
            else {
              m_softGroomed[i] = false;
              m_allSoftgmodes[i] = GROOM_MODE::SD;
            }
          }
        }
      }
    }

    m_groomed = std::valarray<bool>(false,m_nx);
    for(const std::valarray<bool>& cg: m_collGroomed) m_groomed &= cg;
    
    m_L = -std::log(m_logArg);
    m_Lz = -log(m_obss[m_n]->GroomTransitionPoint(p_ampl));
    m_Lsoft = m_L;
    m_Lsoft[m_softGroomed] = m_Lz;
    m_Lsoft[m_groomed and not m_softGroomed] = 0.;
    
    for(const auto& p: m_varParams) {
      msg_Debugging()<<"Resum with variation "<<p.first<<".\n";
      m_params = m_varParams[p.first];
      m_lambda = alphaS()*beta0()*m_L;
      m_lambdaSoft = alphaS()*beta0()*m_Lsoft;
      m_TofLoverA = T(m_lambdaSoft/m_a[0]);

      m_Rp = std::valarray<double>(0., m_nx);

      m_G = Matrix<std::valarray<double>>(4,4,std::valarray<double>(0.,m_nx));
      m_Rexp = Matrix<std::valarray<double>>(4,4,std::valarray<double>(0.,m_nx));
      m_S1 = 0; 
      m_S2 = 0;
      m_P1 = 0;
      m_F2 = 0;
      m_SNGL2 = 0;
      m_EP1 = std::valarray<double> (0., m_nx);
      m_EP2 = std::valarray<double> (0., m_nx);
      m_RAtEnd = std::valarray<double> (0., m_nx);


      
      m_Coll  =  (m_amode&MODE::IGNCOLL)  ? std::valarray<double>(1.,m_nx) : CalcColl(m_Rp,m_G,m_Rexp,m_S1,m_RAtEnd);
      m_Soft  =  (m_amode&MODE::IGNSOFT)  ? std::valarray<double>(1.,m_nx) : CalcS(m_S1,m_S2);
      if(m_amode & MODE::CKINV) {
        double dummy1, dummy2;
        const std::valarray<double> ckinv = std::abs(m_Soft-CalcS(dummy1, dummy2, MODE::CKINV)).apply([](double x)->double {return FuzzyZeroToZero(NaNtoZero(x));});
        const bool failed =  not IsZero(ckinv.max(),1e-6);
        msg_Debugging()<<"Check inversion -> "<<ckinv.max()<<".\n";
        if(failed) {
          msg_Error()<<"Check inversion -> "<<ckinv.max()<<".\n";
          msg_Error()<<"L\t diff\n";
          printLists(msg_Error(),{m_L,ckinv});
        }
      }
      if(m_amode & MODE::CKCOUL) {
        double dummy1, dummy2;
        const std::valarray<double> ckcoul = std::abs(m_Soft-CalcS(dummy1, dummy2, MODE::CKCOUL)).apply([](double x)->double {return FuzzyZeroToZero(NaNtoZero(x));});;
        const bool singletInitial = p_ampl->Leg(0)->Flav().StrongCharge()==0 and p_ampl->Leg(1)->Flav().StrongCharge()==0;
        const bool failed = (singletInitial and not IsZero(ckcoul.max(),1e-6)) or (not singletInitial and IsZero(ckcoul.max(),1e-6));
        msg_Debugging()<<"Check couomb cancellation.\n";
        if(failed) {
          if(singletInitial) {
            msg_Error()<<"Coulomb phases did not cancel despite singlet initial state.\n";
            msg_Error()<<"L\t diff\n";
            printLists(msg_Error(), {m_L,ckcoul});
          }
          else {
            msg_Out()<<"Complete cancellation of coulomb phases with non-singlet II states.\n";
            msg_Out()<<"L\t diff\n";
            printLists(msg_Out(), {m_L,ckcoul});
          }
        }
      }
      m_PDF   =  (m_amode&MODE::IGNPDF)   ? std::valarray<double>(1.,m_nx) : CalcPDF(m_P1);
      m_Ffunc =  (m_amode&MODE::IGNFFUNC) ? std::valarray<double>(1.,m_nx) : CalcF(m_F2);
      m_SNGL  =  (m_amode&MODE::IGNSNGL) ? std::valarray<double>(1.,m_nx) : CalcSNGL(m_SNGL2);
      m_Ep    =  CalcEP(m_EP1, m_EP2);

    

      // m_resNLL[m_n] = ( exp(m_Coll)*m_Soft*m_SNGL*m_PDF*m_Ffunc*m_Ep ).apply([](double x)->double {return std::isnan(x)?0.:x;});
      m_varResNLL[m_n][p.first] = ( exp(m_Coll)*m_Soft*m_SNGL*m_PDF*m_Ffunc*m_Ep ).apply([](double x)->double {return FuzzyZeroToZero(NaNtoZero(x));});

      // some legacy options
      if(!(m_amode & MODE::COLLEXPAND)) {
        msg_Debugging()<<"Ignore coll. expansion.\n";
        m_G = MatrixD(4,4,0);
        m_F2 = 0;
      }
      if(!(m_amode & MODE::SOFTEXPAND)) {
        msg_Debugging()<<"Ignore soft expansion.\n";
        m_S1 = 0;
        m_S2 = 0;
      }
      if(!(m_amode & MODE::PDFEXPAND)) {
        msg_Debugging()<<"Ignore pdf expansion.\n";
        m_P1 = 0;
      }


      std::valarray<double> expCollLO = pow(m_L,2.)*m_G(1,2)  + m_L*m_G(1,1)  + m_G(1,0);
      std::valarray<double> expSoftLO = m_Lsoft * 4./m_a[0] * m_S1;
      std::valarray<double> expPdfLO = m_L * m_P1;
      m_varResExpLO[m_n][p.first] = alphaSBar() * ( expCollLO
                                        + expSoftLO
                                        + expPdfLO
                                        + m_EP1
                                        );
      m_varResExpLO[m_n][p.first] = m_varResExpLO[m_n][p.first].apply([](double x)->double {return FuzzyZeroToZero(NaNtoZero(x));});


  
      m_varResExpNLO[m_n][p.first] = pow(alphaSBar(),2)*pow(m_L,4.) * ( 0.5*pow(m_G(1,2),2.) ); // H24
      m_varResExpNLO[m_n][p.first] += pow(alphaSBar(),2)*pow(m_L,3.) * ( m_G(2,3) + m_G(1,2)*(m_G(1,1) ) ); // H23
      m_varResExpNLO[m_n][p.first] += pow(alphaSBar(),2)*pow(m_L,2.) * ( 0.5*pow(m_G(1,1),2.) + m_G(2,2) +
                                                             m_G(1,2)*m_G(1,0) +
                                                             4.*m_F2*pow(m_Rexp(1,2),2.) ); // H22

      m_varResExpNLO[m_n][p.first] += pow(alphaSBar(),2)*pow(m_L,1.) * ( m_G(2,1) + m_G(1,1)*m_G(1,0) +  4.*m_F2*m_Rexp(1,2)*m_Rexp(1,1) ); // H21
      m_varResExpNLO[m_n][p.first] += pow(alphaSBar(),2)*pow(m_L,0.) * ( m_G(2,0) + 0.5*pow(m_G(1,0),2.) + m_F2*pow(m_Rexp(1,1),2.) ); // H20



      m_varResExpNLO[m_n][p.first] += pow(alphaSBar(),2)*pow(m_Lsoft,2.) * (16./pow(m_a[0],2)*(m_S2 + m_SNGL2) + 4./m_a[0]*m_S1*(2.*M_PI*beta0()/m_a[0]));
      m_varResExpNLO[m_n][p.first] += pow(alphaSBar(),2)*expSoftLO*expCollLO;


      m_varResExpNLO[m_n][p.first] += pow(alphaSBar(),2) * (expCollLO+expSoftLO+expPdfLO)*m_EP1;
      m_varResExpNLO[m_n][p.first] += pow(alphaSBar(),2) * m_EP2;

      m_varResExpNLO[m_n][p.first] = m_varResExpNLO[m_n][p.first].apply([](double x)->double {return FuzzyZeroToZero(NaNtoZero(x));});
    }
    m_params = m_varParams["CENTRAL"];
    // m_resNLL[m_n] = m_varResNLL[m_n]["CENTRAL"];
    // m_resExpLO[m_n] = m_varResExpLO[m_n]["CENTRAL"];
    // m_resExpNLO[m_n] = m_varResExpNLO[m_n]["CENTRAL"];
  }
  CleanUp();
  return 1;
}

void Resum::FinishShowers(double weight, int ntrial) {
  std::map<std::string, double> weights;
  if(m_globaltreat & TREAT::needs_resum) {
    weights = VarWeightsMap();
    weights["CENTRAL"] = weight;
  }
  for(size_t k=0; k<m_obss.size(); k++) {
    // std::map<std::string, double> weights;
    const std::string& tag = m_obss[k]->Tag();
    if(m_resultsNLL[tag] and m_treat[k]&TREAT::needs_resum) {
      // std::map<std::string,std::valarray<double>> NLL;
      // std::map<std::string,std::valarray<double>> expLO;
      // std::map<std::string,std::valarray<double>> expNLO;
      // weights["CENTRAL"] = weight;
      // NLL["CENTRAL"] = m_resNLL[k];
      // expLO["CENTRAL"] = m_resExpLO[k];
      // expNLO["CENTRAL"] = m_resExpNLO[k];
      // for(size_t nv=0; nv<VarNum(); nv++) {
      //   const std::string& vName = VarName(nv);
      //   const double vWeight = VarWeight(nv);
      //   weights[vName] = vWeight;
      //   NLL[vName] = m_varResNLL[k][vName];
      //   expLO[vName] = m_varResExpLO[k][vName];
      //   expNLO[vName] = m_varResExpNLO[k][vName];
      // }
      // m_resultsNLL[tag]->Fill(&weights,&NLL,ntrial);
      // m_resultsNLL[tag]->Fill(m_curChannels);    
      // m_resultsExpLO[tag]->Fill(&weights,&expLO,ntrial);
      // m_resultsExpLO[tag]->Fill(m_curChannels);    
      // m_resultsExpNLO[tag]->Fill(&weights,&expNLO,ntrial);
      // m_resultsExpNLO[tag]->Fill(m_curChannels);    

      m_resultsNLL[tag]->Fill(&weights,&m_varResNLL[k],ntrial);
      m_resultsNLL[tag]->Fill(m_curChannels);    
      m_resultsExpLO[tag]->Fill(&weights,&m_varResExpLO[k],ntrial);
      m_resultsExpLO[tag]->Fill(m_curChannels);    
      m_resultsExpNLO[tag]->Fill(&weights,&m_varResExpNLO[k],ntrial);
      m_resultsExpNLO[tag]->Fill(m_curChannels);
    }
    if(m_resultsLO[tag]) {
      m_resultsLO[tag]->AddZeroPointMCB(ntrial);
      m_resultsLO[tag]->FinishMCB();
    }
    if(m_resultsNLO[tag]) {
      m_resultsNLO[tag]->AddZeroPointMCB(ntrial);
      m_resultsNLO[tag]->FinishMCB();
    }
  }
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
  m_kij.clear();
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

  if(m_obss.size() == 0) return true;


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

  m_globaltreat = TREAT::none;
  for(size_t i=0; i<m_obss.size(); i++) {
    m_resFO[i] = m_obss[i]->Value(ampl);
    m_treat[i] = TREAT::none;
    if(ampl->OrderQCD() == m_obss[i]->ResumQCDorderBorn()) m_treat[i] |= TREAT::born;
    if(ampl->OrderQCD() == m_obss[i]->ResumQCDorderLO())   m_treat[i] |= TREAT::lo;
    if(ampl->OrderQCD() == m_obss[i]->ResumQCDorderNLO())  m_treat[i] |= TREAT::nlo;
    if(m_treat[i]&TREAT::born and m_obss[i]->ResumMult().count(nColoredLegs()) > 0)   
      m_treat[i] |= TREAT::needs_resum;
    m_globaltreat |= m_treat[i];
  }

  if(m_globaltreat & TREAT::needs_resum) {

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

    m_ordered_ids = vector<size_t>(p_ampl->Legs().size());
    for(size_t i=0; i<p_ampl->Legs().size(); i++) {
      m_ordered_ids.at(i) = p_ampl->Leg(cmetric->Map(i))->Id();
    }
  
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

    for(size_t i=0; i<p_ampl->Legs().size()-color_sings; i++) {
      for(size_t j=i+1; j<p_ampl->Legs().size()-color_sings; j++) {
        m_kij.push_back({p_ampl->IdIndex(m_ordered_ids[i]),p_ampl->IdIndex(m_ordered_ids[j])});
      }
    }

  }
  else {
    p_ampl->Delete();
    p_ampl=ampl->Copy();
  }

  Vec4D_Vector moms(p_ampl->Legs().size());
  Flavour_Vector flavs(p_ampl->Legs().size());
  for (size_t i(0);i<p_ampl->Legs().size();++i) {
    moms[i]=i<p_ampl->NIn()?-p_ampl->Leg(i)->Mom():p_ampl->Leg(i)->Mom();
    flavs[i]=i<p_ampl->NIn()?p_ampl->Leg(i)->Flav().Bar():p_ampl->Leg(i)->Flav();
  }
  

  m_curChannels = vector<string>(m_channelAlgs.size());
  for(size_t i=0; i<m_curChannels.size(); i++) {
    m_curChannels[i] = m_channelAlgs[i]->Channel(moms,flavs,2,true);
  }
  
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
                                 std::valarray<double>& xvals)
{
  DEBUG_FUNC(key.Name());
  Observable_Base* obs = RESUM::Observable_Getter::GetObject(key.Name(),key);
  if(obs != nullptr) {
    m_obss.push_back(obs);
    // by convention, xvals are always ordred and there is always an xvalue above or at the max endpoint
    const double ep = obs->MaxEndpoint();
    std::sort(std::begin(xvals),std::end(xvals));
    const size_t nx = xvals.max() < ep ? xvals.size()+1 : xvals.size();
    m_xvals.emplace_back(nx);
    std::valarray<double>& xv = m_xvals.back(); 
    for(size_t i=0; i<nx; i++) {
      if(i<xvals.size()) xv[i] = xvals[i];
      else xv[i] = ep;
    }
    m_resNLL.push_back(valarray<double>(nx));
    m_resExpLO.push_back(valarray<double>(nx));
    m_resExpNLO.push_back(valarray<double>(nx));
    m_varResNLL.push_back(std::map<std::string,valarray<double>>());
    m_varResExpLO.push_back(std::map<std::string,valarray<double>>());
    m_varResExpNLO.push_back(std::map<std::string,valarray<double>>());
    m_ress.push_back({-1,-1});
    msg_Debugging()<<"Added "<<obs->Name()<<" as "<<obs->Tag()<<".\n";
  }
  else  msg_Error()<<"Observable not found: "<<key.Name()<<".\n";
  return obs->Tag();
}

std::string Resum::AddObservable(const RESUM::Observable_Key& key,
                                 std::valarray<double>& xvals,
                                 std::string title, bool calcNLL, 
                                 bool calcLO, bool calcNLO)
{
  DEBUG_FUNC(key.Name());
  Observable_Base* obs = RESUM::Observable_Getter::GetObject(key.Name(),key);
  if(obs != nullptr) {
    calcNLL = calcNLL and HasOrder(obs->ResumQCDorderBorn());
    calcLO = calcLO and HasOrder(obs->ResumQCDorderLO());
    calcNLO = calcNLO and HasOrder(obs->ResumQCDorderNLO());
    m_obss.push_back(obs);
    // by convention, xvals are always ordred and there is always an xvalue above or at the max endpoint
    const double ep = obs->MaxEndpoint();
    std::sort(std::begin(xvals),std::end(xvals));
    const size_t nx = xvals.max() < ep ? xvals.size()+1 : xvals.size();
    m_xvals.emplace_back(nx);
    std::valarray<double>& xv = m_xvals.back(); 
    for(size_t i=0; i<nx; i++) {
      if(i<xvals.size()) xv[i] = xvals[i];
      else xv[i] = ep;
    }
    const std::string& tag = obs->Tag();
    m_treat.push_back(TREAT::none);
    m_resFO.push_back(double());
    m_resultsLO[tag] = calcLO ? std::make_shared<Cumulant>(xvals, title, nullptr, 1, "CENTRAL") : nullptr;
    if(m_resultsLO[tag]) {
      msg_Debugging()<<"Add LO for "<<tag<<".\n";
      m_resultsLO[tag]->InitVarWeights(m_varNames);
      m_resultsLO[tag]->SetVariants(m_channels);
    }
    m_resultsNLO[tag] = calcNLO ? std::make_shared<Cumulant>(xvals, title, nullptr, 1, "CENTRAL") : nullptr;
    if(m_resultsNLO[tag]) {
      msg_Debugging()<<"Add NLO for "<<tag<<".\n";
      m_resultsNLO[tag]->InitVarWeights(m_varNames);
      m_resultsNLO[tag]->SetVariants(m_channels);
    }
    m_resultsNLL[tag] = calcNLL ? std::make_shared<Cumulant>(xvals, title, nullptr, 1, "CENTRAL") : nullptr;
    m_resultsExpLO[tag] = calcNLL ? std::make_shared<Cumulant>(xvals, title, nullptr, 1, "CENTRAL") : nullptr;
    m_resultsExpNLO[tag] = calcNLL ? std::make_shared<Cumulant>(xvals, title, nullptr, 1, "CENTRAL") : nullptr;
    if(m_resultsNLL[tag]) {
      msg_Debugging()<<"Add NLL for "<<tag<<".\n";
      m_resultsNLL[tag]->InitVarWeights(m_varNames);
      m_resultsNLL[tag]->SetVariants(m_channels);
      m_resultsExpLO[tag]->InitVarWeights(m_varNames);
      m_resultsExpLO[tag]->SetVariants(m_channels);
      m_resultsExpNLO[tag]->InitVarWeights(m_varNames);
      m_resultsExpNLO[tag]->SetVariants(m_channels);
      m_resNLL.push_back(valarray<double>(nx));
      m_resExpLO.push_back(valarray<double>(nx));
      m_resExpNLO.push_back(valarray<double>(nx));
      m_varResNLL.push_back(std::map<std::string,valarray<double>>());
      m_varResExpLO.push_back(std::map<std::string,valarray<double>>());
      m_varResExpNLO.push_back(std::map<std::string,valarray<double>>());
    }
    msg_Debugging()<<"Added "<<obs->Name()<<" as "<<obs->Tag()<<".\n";
  }
  else  msg_Error()<<"Observable not found: "<<key.Name()<<".\n";
  return obs->Tag();
}

void Resum::AddChAlg(const std::string& name, const std::vector<std::string> args) {
  m_channelAlgs.emplace_back(ChAlg_Getter::GetObject(name, {name, args}));
  for(const std::string& ch: m_channelAlgs.back()->ChannelNames(true)) {
    m_channels.push_back(ch);
  }
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

double Resum::T(const double& x) const {
  return -1./beta0()/M_PI * std::log(1.-2.*x);
} 

std::valarray<double> Resum::T(const std::valarray<double>& x) const {
  return -1./beta0()/M_PI * std::log(1.-2.*x);
}


std::valarray<double> Resum::CalcS(double& SoftexpNLL_LO, 
                                   double& SoftexpNLL_NLO, 
                                   MODE Check)
{
  //Exception for n_colored = 2
  if(nColoredLegs() == 2) {
    SoftexpNLL_NLO = pow(SoftexpNLL_LO,2.)/2.;
    return std::valarray<double>(1.,m_nx);
  }
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
    for(size_t k=0; k<nColoredLegs(); k++) {
      size_t k_t = (k < 2 ? k : 2*k-1);
      double Cl = flavlabels[k_t]==Flavour(kf_gluon) ? m_params.CA() : m_params.CF();
      Csum += Cl*met;
    }
    msg_Debugging()<<(2.*Tsum+Csum).setFuzzyZeroToZeroInline()<<std::endl;    
  }
    

  //Build Gamma
  MatrixC Gamma(dim, dim, 0);
  for(size_t k=0; k<Tprods.size(); k++) {
    Gamma += m_obss[m_n]->SoftGlobal(p_ampl,m_kij[k].first, m_kij[k].second, s_12)*Tprods[k];
    if(signlabels[2*k]*signlabels[2*k+1] == 1) {
      Gamma -= complex<double>(0,M_PI/2.)*MatrixC(Tprods[k]);
    }
  }

  if(msg_LevelIsDebugging()) {
    msg_Debugging()<<"Gamma = \n"<<Gamma<<"\n";
  }
  if(Check & MODE::CKCOUL) {
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
  SoftexpNLL_LO += 2.*Trace(Hard, Gamma.real())/traceH;
  MatrixC conjGamma = Conjugate(Gamma);
  if(m_mmode & MATCH_MODE::NLO) {
    if((m_amode & MODE::SOFTEXPAND)) {
      SoftexpNLL_NLO *= (SoftexpNLL_LO-SoftexpNLL_NLO/2.);
      SoftexpNLL_NLO += 4.*(Trace(Hard,real(conjGamma*ICmetric*conjGamma))
                            + 2.*Trace(Hard,real(conjGamma*ICmetric*Gamma))
                            + Trace(Hard,real(Gamma*ICmetric*Gamma)))/traceH/8.;
    }
  }
  // Calculate Soft matrix
  
  //  Soft = m_TofLoverA*IGamma;
  std::valarray<MatrixC> eGamma(ICmetric*Gamma.transposeInPlace(),m_nx);
  // for(size_t i=0; i<m_nx; i++) {eGamma[i] *= (m_TofLoverA[i]*eGamma[i]).exponentiateInPlace();}
  for(size_t i=0; i<m_nx; i++) {eGamma[i] *= m_TofLoverA[i];}
  eGamma = exp(eGamma);
  std::valarray<MatrixC> Soft = Conjugate(eGamma)*eGamma;
  
  // if(msg_LevelIsDebugging()) {
  //   msg_Debugging()<<"Soft Matrix = \n"<<Soft.setFuzzyZeroToZero()<<"\n";
  // }
  //Hard-soft contraction 
  std::valarray<double> traceHS = Trace(Soft,Hard*met);
  
 //  if(msg_LevelIsDebugging()) {
 //    msg_Debugging()<<"==============================================" << std::endl;
 //    msg_Debugging()<<"Soft function:" << std::endl;
 //    msg_Debugging()<<"==============================================" << std::endl;   
  //    msg_Debugging()<< "alpha_s: " << as << std::endl;
 //    msg_Debugging()<< "evolution variable t: " << t << std::endl;
 //    msg_Debugging()<< "Log(1/v): " << L << std::endl;
 //    msg_Debugging()<< std::endl;
 //    msg_Debugging()<< "Kinematics" << std::endl;
 //    msg_Debugging()<< *p_ampl << std::endl;
 //    msg_Debugging()<< "Check energy-momentum conservation: ";
 //   Vec4D testEM(0,0,0,0);
 //   for(const auto& l: p_ampl->Legs()) testEM += l->Mom();
 //   msg_Debugging()<< testEM << std::endl;
 //   msg_Debugging()<<"Tr( c H ): " << traceH << std::endl;
 //   msg_Debugging()<<"Softexp lo: " << SoftexpNLL_LO << std::endl;
 //   if((m_amode & MODE::SOFTEXPAND) && (m_mmode & MATCH_MODE::NLO))
 //     msg_Debugging()<<"Softexp nlo: " << SoftexpNLL_NLO << std::endl;
 //   msg_Debugging() <<"Tr( H G Gb ) / Tr( c H ): " << traceHS/traceH << std::endl;
 //   msg_Debugging() << std::endl;
 //   //Print Tprods
 //   size_t k = 0;
 //   for(size_t i = 0; i<nColoredLegs(); i++){
 //     for(size_t j = i+1; j<nColoredLegs(); j++){
 //       const double traceHT = Trace(Tprods[k],Hard);
 //       msg_Debugging() << "T" << i << ".T" << j << "  :  " << std::endl;
 //       msg_Debugging() << "T" << i << " Flavour: " << flavlabels[2*k] << ":  four vec: " << momlabels[2*k] << std::endl;
 //       msg_Debugging() << "T" << j << " Flavour: " << flavlabels[2*k+1] << ":  four vec: " << momlabels[2*k+1] << std::endl;
 //       msg_Debugging() << "log(Qij/Q12): " << log(m_Qij[k]/s_12) << std::endl;
 //       msg_Debugging() << "Qij*Qij: " << std::setprecision(9) << m_Qij[k]*m_Qij[k] << std::endl;
 //       msg_Debugging() <<  "Tr( T.H )/Tr(c.H): " << traceHT/traceH << std::endl;
 //       msg_Debugging() << std::endl;
 //       msg_Debugging()<<Tprods[k]<<std::endl;
 //       msg_Debugging() << std::endl;
 //       k++;
 //     }
 //   }
 //   msg_Debugging()<<"Going to return traceHS/traceH = "<< traceHS/traceH<<std::endl;
 //   msg_Debugging() << std::endl;
 //   msg_Debugging()<<"==============================================" << std::endl;
 //   msg_Debugging()<<"end checks" << std::endl;
 //   msg_Debugging()<<"==============================================" << std::endl;
 //   msg_Debugging()<< std::endl;
 // }
 return traceHS/traceH;
}


double Resum::CalcRpp(const double L, RESUM::GROOM_MODE gmode, double &exp12){
    const double muR2 = p_ampl->MuR2();
    const double beta0 = m_params.beta0(muR2);
    
    double Rpp = 0.;
    
    Poincare cms(p_ampl->Leg(0)->Mom()+p_ampl->Leg(1)->Mom());
    
    Vec4D_Vector moms(p_ampl->Legs().size());
    Flavour_Vector flavs(p_ampl->Legs().size());
    for (size_t i(0);i<p_ampl->Legs().size();++i) {
        moms[i]=i<p_ampl->NIn()?-p_ampl->Leg(i)->Mom():p_ampl->Leg(i)->Mom();
        flavs[i]=i<p_ampl->NIn()?p_ampl->Leg(i)->Flav().Bar():p_ampl->Leg(i)->Flav();
    }
    
    for(size_t i = 2; i<p_ampl->Legs().size(); i++) {
        if(m_deltad[i] == 0) continue;
        msg_Debugging()<<"Calculate radiator for leg "<<i<<".\n";
        double colfac = 0.;
        
        const double as = m_params.alphaS(muR2);
        
        if (p_ampl->Leg(i)->Flav().StrongCharge() == 8) {
            colfac = m_params.CA();
            msg_Debugging()<<"Gluon, Cl = "<<colfac<<".\n";
        }
        else if (abs(p_ampl->Leg(i)->Flav().StrongCharge()) == 3) {
            colfac = m_params.CF();
            msg_Debugging()<<"Quark, Cl = "<<colfac<<".\n";
        } 
        else {
            msg_Debugging()<<"No strong charge, ignoring.\n";
            continue;
        }
        
        
        const double lambda = as*beta0*L;
        
        msg_Debugging()<<"lambda = as*beta0*L = "<<as<<"*"<<beta0<<"*"<<L<<" = "<<lambda<<"\n";
        
        // needed for SD grooming
        const double transp = m_obss[m_n]->GroomTransitionPoint(p_ampl, i);
        msg_Debugging() << "Transition point = " << transp << " zcut = " << m_zcut << "\n";
        
        const double lambdaZ = as*beta0*log(1./transp)/m_a[i];
        const double lambda2 = as*beta0*log(1./2.);
        
        
        if(gmode & GROOM_MODE::SD) {
            Rpp += -colfac*2.*as/M_PI*(2.*lambdaZ+m_beta)/(m_a[i]+m_b[i]-2.*lambda)/(m_a[i]*(1.+m_beta)+m_b[i]-2.*(1.+m_beta)*lambda-2.*m_b[i]*lambdaZ);
            exp12 += -2.*colfac*m_beta/(m_a[i]+m_b[i])/(m_a[i]*(1.+m_beta)+m_b[i]);
        }
        else {
            Rpp += -colfac*2.*as/M_PI/(m_a[i]-2.*lambda)/(m_a[i]+m_b[i]-2.*lambda);
            exp12 += -2./m_a[i] * colfac/(m_a[i]+m_b[i]);
        }
    }    
    return(Rpp);
}

std::valarray<double> Resum::r1(std::valarray<double> lambda, double a, double b) const {
  if(IsZero(b)) return r1_b0(lambda,a);
  return 1./2./M_PI/pow(beta0(),2.)/alphaS()/b*((a-2.*lambda)*log(1.-2.*lambda/a)
                                          -(a+b-2.*lambda)*log(1.-2.*lambda/(a+b)));
}

std::valarray<double> Resum::r1_b0(std::valarray<double> lambda, double a) const {
  return -1./2./M_PI/pow(beta0(),2.)/alphaS()*(2.*lambda/a+log(1.-2.*lambda/a));
}


std::valarray<double> Resum::SD_r1(std::valarray<double> lambda, double lambdaZ, double a, double b) const {
  if(IsZero(b)) return SD_r1_b0(lambda,lambdaZ,a);
  return -1./2./M_PI/pow(beta0(),2)/alphaS()/b * (b/(1.+m_beta) * (1.-2.*lambdaZ)*log(1.-2.*lambdaZ) \
                                              - (a*(1.+m_beta)+b)/(1.+m_beta) * (1.-2*(1.+m_beta)/(a*(1.+m_beta)+b)*lambda - 2.*b/(a*(1.+m_beta)+b)*lambdaZ)*log(1.-2*(1.+m_beta)/(a*(1.+m_beta)+b)*lambda - 2.*b/(a*(1.+m_beta)+b)*lambdaZ) \
                                              + (a+b)*(1.-2./(a+b)*lambda)*log(1.-2./(a+b)*lambda) );
}

std::valarray<double> Resum::SD_r1_b0(std::valarray<double> lambda, double lambdaZ, double a) const {
  return -1./2./M_PI/pow(beta0(),2.)/alphaS()*(m_beta*(2.*lambda/a+log(1.-2.*lambda/a))+2.*lambdaZ+log(1.-2.*lambdaZ)+2.*lambdaZ*(log(1.-2.*lambda/a)-log(1.-2.*lambdaZ)))/(1.+m_beta);
}



std::valarray<double> Resum::r2_cmw(std::valarray<double> lambda, double a, double b) const {
  if(IsZero(b)) return r2_cmw_b0(lambda,a);
  return (K_CMW()/pow(2.*M_PI*beta0(),2.)+log(muR2()/muQ2())/M_PI/beta0()/2.)*((a+b)*log(1.-2.*lambda/(a+b))
                                                               -a*log(1.-2.*lambda/a))/b;
}

std::valarray<double> Resum::r2_beta1(std::valarray<double> lambda, double a, double b) const {
  if(IsZero(b)) return r2_beta1_b0(lambda,a);
  return beta1()/2./M_PI/pow(beta0(),3.)*(a/2.*pow(log(1-2.*lambda/a),2.)
                                   -0.5*(a+b)*pow(log(1.-2.*lambda/(a+b)),2.)
                                   +a*log(1-2.*lambda/a)
                                   -(a+b)*log(1.-2.*lambda/(a+b)))/b;
}

std::valarray<double> Resum::r2_cmw_b0(std::valarray<double> lambda, double a) const {
  return (K_CMW()/pow(2.*M_PI*beta0(),2.)+log(muR2()/muQ2())/M_PI/beta0()/2.)*(log(1.-2.*lambda/a)+2./a*lambda/(1.-2./a*lambda));
}

std::valarray<double> Resum::r2_beta1_b0(std::valarray<double> lambda, double a) const {
  return -beta1()/2./M_PI/pow(beta0(),3.)*(1./2.*pow(log(1-2.*lambda/a),2.)
                                                    +(log(1-2.*lambda/a)+2./a*lambda)/(1.-2.*lambda/a));
}

std::valarray<double> Resum::r1p_b0(std::valarray<double> lambda, double a) const {
  return 2./sqr(a)/(M_PI*beta0())*lambda/(1.-2.*lambda/a);
}


std::valarray<double> Resum::r1p(std::valarray<double> lambda, double a, double b) const {
  if(IsZero(b)) return r1p_b0(lambda,a);
  return 1./b*(T(lambda/a)-T(lambda/(a+b)));
}

std::valarray<double> Resum::SD_r2_cmw(std::valarray<double> lambda, double lambdaZ, double a, double b) const {
  if(IsZero(b)) return SD_r2_cmw_b0(lambda,lambdaZ,a);
  return (K_CMW()/pow(2.*M_PI*beta0(),2.)+log(muR2()/muQ2())/M_PI/beta0()/2.) * (b/(1.+m_beta) * log(1.-2.*lambdaZ) \
                                                                   - (a*(1.+m_beta)+b)/(1.+m_beta) * log(1.-(2.*(1.+m_beta))/(a*(1.+m_beta)+b)*lambda - 2.*b/(a*(1.+m_beta)+b)*lambdaZ) + (a+b)*log(1.-2./(a+b)*lambda))/b;
}

std::valarray<double> Resum::SD_r2_beta1(std::valarray<double> lambda, double lambdaZ, double a, double b) const {
  if(IsZero(b)) return SD_r2_beta1_b0(lambda,lambdaZ,a);
  return -(beta1()/4./M_PI/pow(beta0(),3.)) * (b/(1.+m_beta) * (pow(log(1.-2.*lambdaZ),2)+2.*log(1.-2.*lambdaZ)) \
                                    -(a*(1.+m_beta)+b)/(1.+m_beta) * (pow(log(1.-2.*(1.+m_beta)/(a*(1.+m_beta)+b)*lambda - 2.*b/(a*(1.+m_beta)+b)*lambdaZ),2)+2.*log(1.-2.*(1.+m_beta)/(a*(1.+m_beta)+b)*lambda - 2.*b/(a*(1.+m_beta)+b)*lambdaZ)) + (a+b)*(pow(log(1.-2./(a+b)*lambda),2)+2.*log(1.-2./(a+b)*lambda)) )/b;
}

std::valarray<double> Resum::SD_r1p(std::valarray<double> lambda, double lambdaZ, double a, double b) const {
  if(IsZero(b)) return SD_r1p_b0(lambda,lambdaZ,a);
  return -1./(M_PI*b*beta0())*(log(1.-2.*(1.+m_beta)/(a*(1.+m_beta)+b)*lambda-2.*b/(a*(1.+m_beta)+b)*lambdaZ)-log(1.-2./(a+b)*lambda));
}  
      
std::valarray<double> Resum::SD_r1d(std::valarray<double> lambda, double lambdaZ, double a, double b) const {
  if(IsZero(b)) return SD_r1d_b0(lambda,lambdaZ,a);
  return -1./(M_PI*a*beta0())/(m_beta+1.)*(log(1.-2.*(1.+m_beta)/(a*(1.+m_beta)+b)*lambda-2.*b/(a*(1.+m_beta)+b)*lambdaZ)-log(1.-2.*lambdaZ));
}

std::valarray<double> Resum::SD_r2_cmw_b0(std::valarray<double> lambda, double lambdaZ, double a) const {
  return (K_CMW()/pow(2.*M_PI*beta0(),2.)+log(muR2()/muQ2())/M_PI/beta0()/2.)/(1.+m_beta)*(log(1.-2.*lambdaZ)+m_beta*log(1.-2.*lambda/a)+2.*(m_beta*lambda/a+lambdaZ)/(1.-2.*lambda/a)); 
} 
             
std::valarray<double> Resum::SD_r2_beta1_b0(std::valarray<double> lambda, double lambdaZ, double a) const {
  return -beta1()/4./M_PI/pow(beta0(),3.)*(log(1.-2.*lambdaZ)*(2.+log(1.-2.*lambdaZ))+m_beta*pow(log(1.-2.*lambda/a),2)+2.*(m_beta+2.*lambdaZ)/(1.-2.*lambda/a)*log(1.-2.*lambda/a)+4.*(m_beta*lambda/a+lambdaZ)/(1.-2.*lambda/a))/(1.+m_beta);
} 
             
std::valarray<double> Resum::SD_r1p_b0(std::valarray<double> lambda, double lambdaZ, double a) const {
  return 2./a/(M_PI*beta0())*(m_beta*lambda/a+lambdaZ)/(1.-2.*lambda/a)/(1.+m_beta);
}

std::valarray<double> Resum::SD_r1d_b0(std::valarray<double> lambda, double lambdaZ, double a) const {
  return 1./a/(M_PI*beta0())*(log(1.-2.*lambdaZ)-log(1.-2.*lambda/a))/(1.+m_beta);
}



std::valarray<double> Resum::CalcColl(std::valarray<double>& Rp, Matrix<std::valarray<double>>& G, 
                                      Matrix<std::valarray<double>>& Rexp, double& S1, 
                                      std::valarray<double>& RAtEnd) {

  std::valarray<double> R(0., m_nx);

  for(size_t i = (m_gmode & GROOM_MODE::SD) ? 2 : 0;
      i<p_ampl->Legs().size(); i++) {
    if(m_deltad[i] == 0 or p_ampl->Leg(i)->Flav().StrongCharge() == 0) continue;

    msg_Debugging()<<"Calculate radiator for leg "<<i<<".\n";
    const double colfac = C(i);
    const double hardcoll = B(i);
    const double El = E(i);
    msg_Debugging()<<"Cl = "<<colfac<<", Bl = "<<hardcoll<<", El = "<<El<<".\n";
    const std::valarray<bool>& groomed =  m_collGroomed.at(i);
    const std::valarray<bool>& unGroomed = not groomed;
    const size_t nGroomed = std::count(std::begin(groomed),std::end(groomed),true);
    const size_t nUnGroomed = m_nx - nGroomed;


    // needed for SD grooming
    const double transp = m_obss[m_n]->GroomTransitionPoint(p_ampl, i);
    msg_Debugging() << "Transition point = " << transp << " zcut = " << m_zcut << "\n";
    const double lambdaZ = alphaS()*beta0()*log(1./transp)/m_a[i];
   
    int order = 1;


    if (order>=0) {    
      //LL part
      R[groomed] -= colfac*SD_r1(m_lambda[groomed],lambdaZ,m_a[i],m_b[i]);
      R[unGroomed] -= colfac*r1(m_lambda[unGroomed],m_a[i],m_b[i]);
    } // end LL
    if (order>=1) {
      //NLL part
      std::valarray<double> r2(0.,m_nx);
      std::valarray<double> r1prime(m_nx);
      std::valarray<double> r1dot(m_nx);
      
      r2[groomed] += SD_r2_cmw(m_lambda[groomed],lambdaZ,m_a[i],m_b[i]);
      r2[unGroomed] += r2_cmw(m_lambda[unGroomed],m_a[i],m_b[i]);
      
      r2[groomed] += SD_r2_beta1(m_lambda[groomed],lambdaZ,m_a[i],m_b[i]);
      r2[unGroomed] += r2_beta1(m_lambda[unGroomed],m_a[i],m_b[i]);


      r1prime[groomed] = SD_r1p(m_lambda[groomed],lambdaZ,m_a[i],m_b[i]);
      r1dot[groomed] = SD_r1d(m_lambda[groomed],lambdaZ,m_a[i],m_b[i]);
      
      r1prime[unGroomed] = r1p(m_lambda[unGroomed],m_a[i],m_b[i]);
      
      // add NLL parts to R and Rp
      const double r1p_coeff = m_logdbar[i]-m_b[i]*log(2.0*El/muQ()) + m_logFac;


      R -= colfac*(r2+r1p_coeff*r1prime+hardcoll*T(m_lambda/(m_a[i]+m_b[i])));
      double r1d_coeff = (m_logdbar[i]+m_logFac+m_a[i]*log(2.*El/muQ())-(m_b[i]+(1.+m_beta)*m_a[i])*log(2.*El/s_12)+m_beta*log(2.0*El/muQ()));
      double T_coeff = log(s_12/muQ()); 
      if(!IsZero(m_etamin[i])) {
        r1d_coeff += -(m_b[i]+(1.+m_beta)*m_a[i])*(m_etamin[i]-log(2.*El/s_12));
        T_coeff += -(m_etamin[i]-log(2.*El/s_12));
      }
      r1dot *= colfac*r1d_coeff;
      R[groomed] -= std::valarray<double>(r1dot[groomed]);

      // this is almost but not quite the same being based on lambda_soft, since the collinear transition point might be different
      std::valarray<double> Targ = m_lambda/m_a[i];
      Targ[groomed] = lambdaZ;

      R -= colfac*T_coeff*T(Targ);
      msg_Debugging()<<"r1p_coeff = "<<r1p_coeff<<", T_coeff = "<<T_coeff<<"\n";

      Rp += r1prime*colfac;        
    }

    msg_Debugging()<<"Now calculate expansion.\n";

    G(1,2)[groomed] += std::valarray<double>(-2.*colfac*m_beta/(m_a[i]+m_b[i])/(m_a[i]*(1.+m_beta)+m_b[i]),nGroomed);
    G(1,1)[groomed] += std::valarray<double>( -4.*colfac*(log(1./transp)/m_a[i]/(m_a[i]*(1.+m_beta)+m_b[i]) + hardcoll/(m_a[i]+m_b[i]) + m_beta/(m_a[i]*(1.+m_beta)+m_b[i])/(m_a[i]+m_b[i])*(m_logdbar[i]-m_b[i]*log(2.0*El/muQ())+m_logFac)+(m_logdbar[i]+m_logFac+m_a[i]*log(2.*El/muQ())-(m_b[i]+(1.+m_beta)*m_a[i])*log(2.*El/s_12)+m_beta*log(2.0*El/muQ()))/m_a[i]/(m_a[i]*(1.+m_beta)+m_b[i])), nGroomed);
    G(1,0)[groomed] += std::valarray<double>( 4.*colfac*(sqr(log(1./transp))/2.0/m_a[i]/(m_a[i]*(1.+m_beta)+m_b[i]) - (1./m_a[i]/(m_a[i]*(1.+m_beta)+m_b[i])*(m_logdbar[i]-m_b[i]*log(2.0*El/muQ())+m_logFac)-(m_logdbar[i]+m_logFac+m_a[i]*log(2.*El/muQ())-(m_b[i]+(1.+m_beta)*m_a[i])*log(2.*El/s_12)+m_beta*log(2.0*El/muQ()))/m_a[i]/(m_a[i]*(1.+m_beta)+m_b[i]))*log(1./transp)), nGroomed );
        
    Rexp(1,2)[groomed] += std::valarray<double>( -2.*colfac*m_beta/(m_a[i]+m_b[i])/(m_a[i]*(1.+m_beta)+m_b[i]), nGroomed );
    Rexp(1,1)[groomed] += std::valarray<double>( -4.*colfac*log(1./transp)/m_a[i]/(m_a[i]*(1.+m_beta)+m_b[i]), nGroomed );

               
    G(2,3)[groomed] += std::valarray<double>(-8.*M_PI*beta0()*colfac/3.*m_beta*(2.*(m_beta+1.)*m_a[i]+(m_beta+2.)*m_b[i])/sqr((m_a[i]+m_b[i])*(m_a[i]*(1.+m_beta)+m_b[i])), nGroomed);
    G(2,2)[groomed] += std::valarray<double>(-8.*M_PI*beta0()*colfac*((m_beta+1.)/m_a[i]/sqr(m_a[i]*(1.+m_beta)+m_b[i])*log(1./transp)
                                          +m_beta*(K_CMW()/2./M_PI/beta0()+log(muR2()/muQ2()))/2./(m_a[i]+m_b[i])/(m_a[i]*(1.+m_beta)+m_b[i])
                                          +m_beta*(2.*(m_beta+1.)*m_a[i]+(m_beta+2.)*m_b[i])/sqr((m_a[i]+m_b[i])*(m_a[i]*(1.+m_beta)+m_b[i]))*(m_logdbar[i]-m_b[i]*log(2.*El/muQ())+m_logFac)
                                                                      +(1.+m_beta)/m_a[i]/sqr(m_a[i]*(1.+m_beta)+m_b[i])*(m_logdbar[i]+m_logFac+m_a[i]*log(2.*El/muQ())-(m_b[i]+(1.+m_beta)*m_a[i])*log(2.*El/s_12)+m_beta*log(2.0*El/muQ()))
                                                                      +hardcoll/sqr(m_a[i]+m_b[i])), nGroomed );
    G(2,1)[groomed] += std::valarray<double>(-8.*M_PI*beta0()*colfac*(m_b[i]/sqr(m_a[i]*(m_a[i]*(1.+m_beta)+m_b[i]))*sqr(log(1./transp)) 
                                          +( (K_CMW()/2./M_PI/beta0()+log(muR2()/muQ2()))/m_a[i]/(m_a[i]*(1.+m_beta)+m_b[i])
                                           +2.*(m_beta+1.)/m_a[i]/sqr(m_a[i]*(1.+m_beta)+m_b[i])*(m_logdbar[i]-m_b[i]*log(2.*El/muQ())+m_logFac)
                                             +2.*m_b[i]/sqr(m_a[i]*(m_a[i]*(1.+m_beta)+m_b[i]))* (m_logdbar[i]+m_logFac+m_a[i]*log(2.*El/muQ())-(m_b[i]+(1.+m_beta)*m_a[i])*log(2.*El/s_12)+m_beta*log(2.0*El/muQ())))*log(1./transp)), nGroomed );
    G(2,0)[groomed] += std::valarray<double>(-8.*M_PI*beta0()*colfac*(-(m_a[i]*(1.+m_beta)+2.*m_b[i])/3./sqr(m_a[i]*(m_a[i]*(1.+m_beta)+m_b[i]))*pow(log(1./transp),3)
                                          +( -(K_CMW()/2./M_PI/beta0()+log(muR2()/muQ2()))/2./m_a[i]/(m_a[i]*(1.+m_beta)+m_b[i])
                                           +m_b[i]/sqr(m_a[i]*(m_a[i]*(1.+m_beta)+m_b[i]))*(m_logdbar[i]-m_b[i]*log(2.*El/muQ())+m_logFac)
                                             -(m_a[i]*(1.+m_beta)+2.*m_b[i])/sqr(m_a[i]*(m_a[i]*(1.+m_beta)+m_b[i]))*(m_logdbar[i]+m_logFac+m_a[i]*log(2.*El/muQ())-(m_b[i]+(1.+m_beta)*m_a[i])*log(2.*El/s_12)+m_beta*log(2.0*El/muQ())))*sqr(log(1./transp))), nGroomed );

    if(!IsZero(m_etamin[i])) {
      G(1,1)[groomed] += std::valarray<double>( 4.*colfac*(m_b[i]+(1.+m_beta)*m_a[i])*(m_etamin[i]-log(2.*El/s_12))/m_a[i]/(m_a[i]*(1.+m_beta)+m_b[i]), nGroomed );
      G(1,0)[groomed] += std::valarray<double>( -4.*colfac*(m_b[i]+(1.+m_beta)*m_a[i])*(m_etamin[i]-log(2.*El/s_12))/m_a[i]/(m_a[i]*(1.+m_beta)+m_b[i])*log(1./transp), nGroomed );
      G(2,2)[groomed] += std::valarray<double>( 8.*M_PI*beta0()*colfac*(1.+m_beta)/m_a[i]/sqr(m_a[i]*(1.+m_beta)+m_b[i])*(m_b[i]+(1.+m_beta)*m_a[i])*(m_etamin[i]-log(2.*El/s_12)), nGroomed );
      G(2,1)[groomed] += std::valarray<double>( 8.*M_PI*beta0()*colfac*2.*m_b[i]/sqr(m_a[i]*(m_a[i]*(1.+m_beta)+m_b[i]))*(m_b[i]+(1.+m_beta)*m_a[i])*(m_etamin[i]-log(2.*El/s_12))*log(1./transp), nGroomed );
      G(2,0)[groomed] += std::valarray<double>( -8.*M_PI*beta0()*colfac*(m_a[i]*(1.+m_beta)+2.*m_b[i])/sqr(m_a[i]*(m_a[i]*(1.+m_beta)+m_b[i]))*(m_b[i]+(1.+m_beta)*m_a[i])*(m_etamin[i]-log(2.*El/s_12))*sqr(log(1./transp)), nGroomed );
          
      G(1,0)[groomed] += std::valarray<double>( 4.*colfac*(m_etamin[i]-log(2.*El/s_12))/m_a[i] * log(1./transp), nGroomed );
      G(2,0)[groomed] += std::valarray<double>( 8.*M_PI*beta0()*colfac*(m_etamin[i]-log(2.*El/s_12))/sqr(m_a[i]) * sqr(log(1./transp)), nGroomed );
    }

 
    RAtEnd[groomed] += std::valarray<double>( colfac*log(1.-2.*m_b[i]*lambdaZ/(m_a[i]*(1.+m_beta)+m_b[i]))/M_PI/m_b[i]/beta0(), nGroomed );

    const double r2_beta1AtEnd = -alphaS()*beta1()/M_PI/beta0()/beta0()*m_a[i]*(log(1.-2.*m_b[i]*lambdaZ/(m_a[i]*(1.+m_beta)+m_b[i]))+2.*m_b[i]*lambdaZ/(m_a[i]*(1.+m_beta)+m_b[i]))/(1.-2.*m_b[i]*lambdaZ/(m_a[i]*(1.+m_beta)+m_b[i]));
    const double r2_cmwAtEnd = 2.*alphaS()*beta0()*(K_CMW()/pow(2.*M_PI*beta0(),2.)+log(muR2()/muQ2())/M_PI/beta0()/2.)*(1./(1.-2.*m_b[i]*lambdaZ/(m_a[i]*(1.+m_beta)+m_b[i]))-1.);
    const double r2_hardcollAtEnd = 2.*alphaS()/M_PI/(m_a[i]+m_b[i]);
    const double r1pAtEnd = -2.*alphaS()/M_PI/m_b[i]*(1./(m_a[i]+m_b[i])-(1.+m_beta)/(m_a[i]*(1.+m_beta)+m_b[i])/(1.-2.*m_b[i]*lambdaZ/(m_a[i]*(1.+m_beta)+m_b[i])));
    const double r1dAtEnd = 2.*alphaS()/M_PI/m_a[i]/(m_a[i]*(1.+m_beta)+m_b[i])/(1.-2.*m_b[i]*lambdaZ/(m_a[i]*(1.+m_beta)+m_b[i]));
        
    const double r2AtEnd=1./m_b[i]*(r2_cmwAtEnd+r2_beta1AtEnd)+m_logFac*r1pAtEnd;
        


    RAtEnd[groomed] += std::valarray<double>( -colfac*(r2AtEnd+r1pAtEnd*(m_logdbar[i]-m_b[i]*log(2.0*El/muQ())) +r1dAtEnd*(m_logdbar[i]+m_logFac+m_a[i]*log(2.*El/muQ())-(m_b[i]+(1.+m_beta)*m_a[i])*log(2.*El/s_12)+m_beta*log(2.0*El/muQ())) +hardcoll*r2_hardcollAtEnd), nGroomed );        

    msg_Debugging()<<"Done expansion for grooming.\n";


    RAtEnd[unGroomed] += std::valarray<double>(-2./M_PI*alphaS()*(colfac) * (hardcoll/(m_a[i]+m_b[i]) + 1./m_a[i]/(m_a[i]+m_b[i])*(m_logdbar[i]-m_b[i]*log(2.0*El/muQ())+m_logFac)), nUnGroomed);
        
        
    G(1,2)[unGroomed] += std::valarray<double>( -2./m_a[i] * colfac/(m_a[i]+m_b[i]), nUnGroomed );
    G(1,1)[unGroomed] += std::valarray<double>( -colfac*(4.*hardcoll/(m_a[i]+m_b[i]) + 4./(m_a[i]*(m_a[i]+m_b[i]))*(m_logdbar[i]-m_b[i]*log(2.*El/muQ())+m_logFac)), nUnGroomed );

    Rexp(1,2)[unGroomed] += std::valarray<double>( -2./m_a[i] * colfac/(m_a[i]+m_b[i]), nUnGroomed );
        
    G(2,3)[unGroomed] += std::valarray<double>( -8.*M_PI*beta0()/3./pow(m_a[i],2) * colfac * (2.*m_a[i]+m_b[i])/pow(m_a[i]+m_b[i],2), nUnGroomed );
    G(2,2)[unGroomed] += std::valarray<double>( -colfac*(8.*M_PI*beta0() * (hardcoll/pow(m_a[i]+m_b[i],2) + 
                                                             (2.*m_a[i]+m_b[i])/pow(m_a[i]*(m_a[i]+m_b[i]),2)*(m_logdbar[i]-m_b[i]*log(2.*El/muQ())+m_logFac)) \
                                          +2.*(K_CMW()+M_PI*beta0()*2.*log(muR2()/muQ2()))/m_a[i]/(m_a[i]+m_b[i])), nUnGroomed );

    if(!IsZero(m_etamin[i])) {
      RAtEnd[unGroomed] += std::valarray<double>( 2./M_PI*alphaS()*colfac*(m_etamin[i]-log(2.*El/s_12))/m_a[i], nUnGroomed );
            
      G(1,1)[unGroomed] += std::valarray<double>( 4.*colfac*(m_etamin[i]-log(2.*El/s_12))/m_a[i], nUnGroomed );
      G(2,2)[unGroomed] += std::valarray<double>( 8.*M_PI*beta0()*colfac*(m_etamin[i]-log(2.*El/s_12))/sqr(m_a[i]), nUnGroomed );
    } 
    msg_Debugging()<<"Done expansion without grooming.\n";

    S1 += -colfac*log(s_12/muQ());
  }
  msg_Debugging()<<"Done looping over legs.\n";
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



std::valarray<double> Resum::CalcPDF(double &PDFexp) {
  if(m_gmode & GROOM_MODE::SD) {
    PDFexp = 0;
    msg_Debugging()<<"Ignoring pdf function for groomed observables.\n";
    return std::valarray<double>(1.,m_nx);
  }
  msg_Debugging()<<"Calculate pdf contribution, no grooming assumed.\n";

  double old_pdffac = 1.;
  std::valarray<double> new_pdffac(1., m_nx);

  const double scale= muF2();
  msg_Debugging() << "scale before: " << scale << "\n";

  for (size_t i=0; i<2; i++) {
    if(m_deltad[i] == 0 or p_ampl->Leg(i)->Flav().IsLepton()) continue;

    const double x = i==0 ? (-p_ampl->Leg(i)->Mom()).PPlus()/rpa->gen.PBeam(0).PPlus() : (-p_ampl->Leg(i)->Mom().PMinus())/rpa->gen.PBeam(1).PMinus();

    //original PDF
    m_params.PDF(i)->Calculate(x,scale);

    const double z = x+(1.0-x)*m_rn[i];

    msg_Debugging()<<"Calculate PDF expansion with z = "<<z<<".\n";
    PDFexp += -2.0/(m_a[i]+m_b[i])*CollinearCounterTerms(i,p_ampl->Leg(i)->Flav().Bar(),-p_ampl->Leg(i)->Mom(),z,scale);

    const double fb = m_params.PDF(i)->GetXPDF(p_ampl->Leg(i)->Flav().Bar());
    old_pdffac *= fb;

    //new PDF scale
    const std::valarray<double> newscale = pow(exp(-m_L),2./(m_a[i]+m_b[i]))*scale;
    for(size_t j=0; j<m_nx; j++) {
      if (newscale[j]<m_params.PDF(i)->Q2Min()) {
        //freeze PDF at Q2Min 
        m_params.PDF(i)->Calculate(x,m_params.PDF(i)->Q2Min());
      }
      else {
        m_params.PDF(i)->Calculate(x,newscale[j]);
      }
      new_pdffac[j] *= m_params.PDF(i)->GetXPDF(p_ampl->Leg(i)->Flav().Bar());      
    }
  }
  return new_pdffac/old_pdffac;
}


std::valarray<double> Resum::CalcF(double& FexpNLL_NLO) {
  std::valarray<double> ret(m_nx);
  for(size_t i=0; i<m_nx; i++) {
    ret[i] = !std::isnan(m_Rp[i]) ? m_F(m_Rp[i],FexpNLL_NLO) : m_F(0.,FexpNLL_NLO);
  }
  return ret;
}

std::valarray<double> Resum::CalcSNGL(double& SnglExpNLL_NLO) {
  std::valarray<double> ret(m_nx);
  for(size_t i=0; i<m_nx; i++)  {
    ret[i] = !std::isnan(m_TofLoverA[i]) ?  m_Sngl(m_TofLoverA[i],SnglExpNLL_NLO) : m_Sngl(0.,SnglExpNLL_NLO);
  }
  return ret;
}

std::valarray<double> Resum::CalcEP( std::valarray<double>&  EPexpNLL_LO,  std::valarray<double>& EPexpNLL_NLO) {
  std::valarray<double> ret(1.,m_nx);
  size_t idx = m_nx-1;
  if(m_mmode & MATCH_MODE::ADD) {
    if(m_gmode & GROOM_MODE::SD) {
      THROW(fatal_error,"Matching mode 'ADD' not valid for grooming.");
    }
    ret *= exp(-alphaSBar()*m_epRatio*(4./m_a[0]*m_S1*m_Lsoft + m_L*(m_P1 + m_G(1,1)[idx])));
    EPexpNLL_LO -=  m_epRatio*(4./m_a[0]*m_S1*m_Lsoft + m_L*(m_P1 + m_G(1,1)[idx]));
    EPexpNLL_NLO += pow(m_epRatio*(4./m_a[0]*m_S1*m_Lsoft + m_L*(m_P1 + m_G(1,1)[idx])),2)/2.; 
  }
  if(m_mmode & MATCH_MODE::DERIV) {
    ret *= exp(-m_Coll[idx])/m_Ffunc[idx];
    EPexpNLL_LO -= std::valarray<double>(m_G(1,0)[idx],m_nx);
    EPexpNLL_NLO -= std::valarray<double>(m_G(2,0)[idx] - 0.5*pow(m_G(1,0)[idx],2) + m_F2*pow(m_Rexp(1,1)[idx],2),m_nx);

    const std::valarray<double> coeff = exp(-m_epRatio*(alphaSBar()*(4./m_a[0]*m_S1*m_Lsoft + m_L*m_P1) + m_L*m_RAtEnd[idx]));
    ret *= coeff;
    EPexpNLL_LO -=  m_epRatio*(4./m_a[0]*m_S1*m_Lsoft + m_L*(m_P1 + m_G(1,1)[idx]));
    EPexpNLL_NLO -= m_epRatio*m_L*m_G(2,1)[idx];

    EPexpNLL_NLO +=  pow(m_epRatio*(4./m_a[0]*m_S1*m_Lsoft + m_L*(m_P1 + m_G(1,1)[idx])),2)/2.;
    EPexpNLL_NLO += m_epRatio*(4./m_a[0]*m_S1*m_Lsoft + m_L*(m_P1 + m_G(1,1)[idx]))*m_G(1,0)[idx];

    if(m_gmode & GROOM_MODE::SD) {
      double temp;
      const double Rpp0 = CalcRpp(0., GROOM_MODE::SD, temp);
      const double Fp0 = (GAMMA_E+DiGamma(1.+m_Rp[idx]))*Rpp0;
      const std::valarray<double> fep = exp(-m_epRatio* m_L * Fp0);
      ret *= fep;
    }
    EPexpNLL_NLO -= m_epRatio*m_L * 4.*m_F2*m_Rexp(1,2)[idx]*m_Rexp(1,1)[idx];
  }
  return ret;
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
