#include "RRatios/RRatios.H"

#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Data_Reader.H"
#include <set>

#include "YODA/Scatter2D.h"
#include "YODA/WriterYODA.h" 


using namespace RESUM;
using namespace PDF;
using namespace ATOOLS;
using namespace MODEL;
using PHASIC::Process_Base;
using std::vector;
using std::string;


inline bool StringToBool(const string& val,
                         const std::set<string>& False = {"0", "False", "false",
                                                            "FALSE", "Off", "off",
                                                            "OFF", "No", "no", "NO"}) {
  return False.find(val)==False.end();
}

RRatios::RRatios(ISR_Handler *const /*isr*/,
                 Model_Base *const model):
  Shower_Base("RRatios")
{
  p_as=(Running_AlphaS*)model->GetScalarFunction("alpha_S");
  
  Data_Reader read(" ",";","#","=");
  m_amode=read.GetValue<int>("RESUM_MODE",0);
  rpa->gen.SetVariable("SCALES", read.GetValue<string>("SCALES", "VAR{sqr(91.188)}"));
  rpa->gen.SetVariable("RESUM::pre_calc", read.GetValue<string>("RESUM::pre_calc", "pre_calc"));
  rpa->gen.SetVariable("RESUM::RRatio_Mode", read.GetValue<string>("RESUM::RRatio_Mode", "Random"));
  const string& energy = read.GetValue<string>("RESUM::RRatio_Energy", "None");
  if(energy != "None") {
    rpa->gen.SetVariable("RESUM::RRatio_Mode", "Starlike");
    if(energy != "FromAmplitude") {
      m_E = read.GetValue<double>("RESUM::RRatio_Energy",-1);
    }
  }
  rpa->gen.SetVariable("RESUM::RRatio_Ratio_name",
                       read.GetValue<string>("RESUM::RRatio_Ratio_name", "Ratio"));
  rpa->gen.SetVariable("RESUM::RRatio_MEn_name",
                       read.GetValue<string>("RESUM::RRatio_MEn_name", "MEn"));
  rpa->gen.SetVariable("RESUM::RRatio_MEnp1_name",
                       read.GetValue<string>("RESUM::RRatio_MEnp1_name", "MEnp1"));
  m_plotRatio = StringToBool(read.GetValue<string>("RESUM::RRatio_plotRatio", "True"));
  m_plotMEn = StringToBool(read.GetValue<string>("RESUM::RRatio_plotMEn", "True"));
  m_plotMEnp1 = StringToBool(read.GetValue<string>("RESUM::RRatio_plotMEnp1", "True"));
  
  m_lambda = read.GetValue<double>("RESUM::RRatio_lambda", 0.95);
  m_cutoff = read.GetValue<double>("RESUM::RRatio_cutoff", 1e-3);
  if (rpa->gen.Variable("SHOWER_GENERATOR")=="") {
    rpa->gen.SetVariable("SHOWER_GENERATOR",ToString(this));
  }
}

RRatios::~RRatios()
{
  while (m_cmetrics.size()>0) {
    delete m_cmetrics.begin()->second;
    m_cmetrics.erase(m_cmetrics.begin());
  }
}

inline kf_code new_flavour(kf_code fl1, kf_code fl2) {
  // TODO: is this correct in all cases?
  if(fl1==21) return fl2;
  if(fl2==21) return fl1;
  return 21;
}

int RRatios::PerformShowers()
{
  DEBUG_FUNC(this);
  // first check we got everything we need
  if(p_ampl_np1==nullptr) THROW(fatal_error,"No process info for n+1.");
  if(p_ampl_n==nullptr) THROW(fatal_error,"No process info for n.");

  if(p_emit==nullptr) THROW(fatal_error,"Emitter not set.");
  if(p_spect==nullptr) THROW(fatal_error,"Spectator not set.");
  if(p_soft==nullptr) THROW(fatal_error,"Soft not set.");

  if(p_cmetric_np1==nullptr) THROW(fatal_error,"No metric for n+1.");
  if(p_cmetric_n==nullptr) THROW(fatal_error,"No metric for n.");


  msg_Debugging()<<*p_ampl_n<<"\n";
  msg_Debugging()<<*p_ampl_np1<<"\n";  
  
  // now rescale one of the momenta
  const Vec4D& emit = p_emit_n->Mom();
  const Vec4D& spect = p_spect_n->Mom();

  const MatrixD& metric_np1 = p_cmetric_np1->CMetric();
  const MatrixD& metric_n = p_cmetric_n->CMetric();


  const MatrixD& pref_np1 = p_cmetric_np1->PrefMatrix();
  const MatrixD& pref_n = p_cmetric_n->PrefMatrix(); 
  
  const MatrixD& trafo_np1 = p_cmetric_np1->TransformationMatrix();
  const MatrixD& trafo_n = p_cmetric_n->TransformationMatrix();

  const MatrixD& trafo_np1_T = Transpose(trafo_np1);
  const MatrixD& trafo_n_T = Transpose(trafo_n);

  const size_t dim_np1 = trafo_np1.numCols();  
  const size_t dim_n = trafo_n.numCols();

  std::vector<MatrixD> Tprods(p_cmetric_n->Tprods().size());
  for(size_t i=0; i<p_cmetric_n->Tprods().size(); i++) {
    Tprods[i] = p_cmetric_n->Tprods().at(i);
  }

  YODA::Scatter2D plot("/ratio/ratio","/ratio/ratio");
  
  YODA::Scatter2D plot_n("/ME/ME","/ME/ME");
  YODA::Scatter2D plot_np1("/ME/ME","/ME/ME");

  for(double cut=m_lambda; cut>m_cutoff; cut*=m_lambda) {
    const Vec4D& soft = m_lambda*p_soft->Mom();
    const double eps = emit*soft/(spect*(emit-soft));
    const Vec4D& p1 = emit-soft + eps * spect;
    const Vec4D& p3 = (1.-eps)*spect;


    msg_Debugging()<<"New momenta:\n"<<p1<<"\n"<<soft<<"\n"<<p3<<"\n\n";
  
    p_emit->SetMom(p1);
    p_soft->SetMom(soft);
    p_spect->SetMom(p3);

    
    MatrixD H_np1 = MatrixC(m_comix.ComputeHardMatrix(p_ampl_np1,
                                                      p_cmetric_np1->Perms()),
                            dim_np1, dim_np1, 0).real();
    H_np1.data() *= pref_np1.data();
    if(p_cmetric_np1->hasTrafo()) H_np1 = trafo_np1*H_np1*trafo_np1_T;
    
    m_comix.Reset();
    msg_Debugging()<<"Read in matrix for n+1 process.\n";
    msg_Debugging()<<H_np1<<std::endl;

    MatrixD H_n = MatrixC(m_comix.ComputeHardMatrix(p_ampl_n,p_cmetric_n->Perms()),
                          dim_n, dim_n, 0).real();
    H_n.data() *= pref_n.data();
    if(p_cmetric_n->hasTrafo()) H_n = trafo_n*H_n*trafo_n_T;
    m_comix.Reset(); //TODO: do I need these resets?
    msg_Debugging()<<"Read in matrix for n process.\n";
    msg_Debugging()<<H_n<<std::endl;

    msg_Debugging()<<"Amplitude for n+1: "<<*p_ampl_np1<<"\n";
    msg_Debugging()<<"Amplitude for n: "<<*p_ampl_n<<"\n";

    msg_Debugging() << "Tr(c_n+1 * H_n+1) = "<< Trace(metric_np1,H_np1)<<"\n";
    msg_Debugging() << "Tr(c_n * H_n) = "<< Trace(metric_n,H_n)<<"\n";

    

    MatrixD Gamma(Tprods.at(0).numRows(), Tprods.at(0).numCols());

    // Do all insertions between -coloured- legs
    size_t i=0;
    for(size_t t=0; t<p_ampl_n->Legs().size(); t++){

      // Check if leg t is coloured
      Flavour flav_t = p_ampl_n->IdLeg(m_ordered_ids.at(t))->Flav();
      if (flav_t.IsGluon() || flav_t.IsQuark()) {
	for(size_t r=t+1; r<p_ampl_n->Legs().size(); r++) {

	  // Check if leg r is coloured
	  Flavour flav_r = p_ampl_n->IdLeg(m_ordered_ids.at(r))->Flav();
	  if (flav_r.IsGluon() || flav_r.IsQuark()) {
	    msg_Debugging()<<"Insertion between legs "
                       <<t<<" -> "<<ID(m_ordered_ids.at(t))<<", "
                       <<r<<" -> "<< ID(m_ordered_ids.at(r))<<".\n";
	    Vec4D pt = p_ampl_n->IdLeg(m_ordered_ids.at(t))->Mom();
	    Vec4D pr = p_ampl_n->IdLeg(m_ordered_ids.at(r))->Mom();
	    // incoming legs have negative momentum
	    // TODO: is there a more natural way to detect them?
	    if(pt[0]<0) pt *= -1;
	    if(pr[0]<0) pr *= -1;
	    const double eikonal = pt*pr / ((pt*soft)*(pr*soft));
	    msg_Debugging()<<"T-product: \n"<<Tprods[i]<<"\n";
	    msg_Debugging()<<"p_t = "<<pt<<"\n";
	    msg_Debugging()<<"p_r = "<<pr<<"\n";
	    msg_Debugging()<<"p_soft = "<<soft<<"\n";
	    msg_Debugging()<<"eikonal = "<<eikonal<<"\n\n";
	    // factor 2 due to symmetrisation of sum
	    Gamma -= 2.*eikonal*Tprods.at(i);
	    i++;
	  }
	}
      }
    }
    /* msg_Debugging()<<"Gamma:\n"<<Gamma<<"\n"; */

    // TODO: is that always the scale we want?
    const double g2 =  4.*M_PI*p_as->AlphaS(p_ampl_np1->MuR2());
    const double TrcH_np1 = Trace(metric_np1,H_np1);

    // we actually want H*metric*colour-change-matrices, but our Tproducts are
    // inverse_metric*colour-change-matrix, so this is correct
    const double TrHG = Trace(H_n,Gamma);
    const double ratio = (g2*TrHG) / TrcH_np1;
    msg_Debugging()<<"softness: "<<cut<<" -> ratio = "<<ratio<<"\n";
    plot.addPoint(cut, ratio);
    plot_n.addPoint(cut,g2*TrHG);
    plot_np1.addPoint(cut,TrcH_np1);
  }
  msg_Debugging()<<"Writing plots to file...";
  if(m_plotRatio) {
    YODA::WriterYODA::write(rpa->gen.Variable("RESUM::RRatio_Ratio_name")
                            +"_"+std::to_string(m_count)+".yoda" ,plot);
  }
  if(m_plotMEn) {
    YODA::WriterYODA::write(rpa->gen.Variable("RESUM::RRatio_MEn_name")
                            +"_"+std::to_string(m_count)+".yoda" ,plot_n);
  }
  if(m_plotMEnp1) {
    YODA::WriterYODA::write(rpa->gen.Variable("RESUM::RRatio_MEnp1_name")
                            +"_"+std::to_string(m_count)+".yoda" ,plot_np1);
  }
  msg_Debugging()<<"done!\n";
  m_count++;
  CleanUp();
  return 1;
}


int RRatios::PerformDecayShowers()
{
  DEBUG_FUNC(this);
  return PerformShowers();
}

bool RRatios::ExtractPartons(Blob_List *const bl)
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

void RRatios::CleanUp()
{
  m_comix.Reset();
}

Cluster_Definitions_Base *RRatios::GetClusterDefinitions()
{
  return nullptr;
}

void RRatios::MinimallyCollinearFinalState(Cluster_Amplitude* ampl) {
  int n = ampl->Legs().size();
  double E = m_E>0 ? m_E:(ampl->Leg(0)->Mom()+ampl->Leg(1)->Mom()).Abs()/2.;
  
  vector<Vec4D> momenta = {{E,0.,0.,E}, {E,0.,0.,-E}};  
  for (int i=0; i<n-2; ++i) {
    momenta.push_back({2*E/(n-2),
          2*E/(n-2)*cos(M_PI*(1./8.+i*2./(n-2))),
          2*E/(n-2)*sin(M_PI*(1./8.+i*2./(n-2))),
          0.});
  }
  SetMomenta(ampl,momenta);
  msg_Debugging()<<"Amplitude after setting into minimally collinear final state:\n";
  msg_Debugging()<<*ampl<<"\n";
}

bool RRatios::PrepareShower
(Cluster_Amplitude* ampl,const bool & /*soft*/)
{
  DEBUG_FUNC(this);
  DEBUG_VAR(ampl->Proc<Process_Base>());

  p_ampl_np1=ampl->Copy();

  // if we want we can reset momenta here
  if(rpa->gen.Variable("RESUM::RRatio_Mode") == "Starlike") {
    MinimallyCollinearFinalState(p_ampl_np1);
  }

  ATOOLS::Cluster_Amplitude* tmp=p_ampl_np1->Copy();
  tmp->SetNIn(0);  
  string pname=Process_Base::GenerateName(tmp);

  
  Process_Base::SortFlavours(tmp);

  CMetric_Base* cmetric;

  n_g=0;
  n_q=0;
  n_aq=0;
   
  for (size_t i=0;i<tmp->Legs().size();i++) {
    Flavour flav = tmp->Leg(i)->Flav();
    if (i<2) flav=flav.Bar();
    if (flav==Flavour(kf_gluon)) n_g++;
    if (flav.IsQuark() && !flav.IsAnti()) n_q++;
    if (flav.IsQuark() && flav.IsAnti()) n_aq++;
  }

  //number of color singlets
  color_sings = tmp->Legs().size() - (n_g + n_aq + n_q);
   
  string name= ToString(n_g)+"_G_"+ToString(n_q)+"_Q_"+ToString(n_aq)+"_AQ";
  
  msg_Debugging()<<" Found process with "<<n_g<<" Gluons, "<<n_q<<" Quarks and "<<n_aq<<" Anti-Quarks : "<<name<<"\n";

  //search for metric in m_cmetrics, 
  //otherwise get it from available getters
  CMetric_Map::const_iterator CMiter=m_cmetrics.find(pname);
  //get_Cbasis(name)
  if (CMiter!=m_cmetrics.end()) {
    cmetric=CMiter->second;
    msg_Debugging()<<" found metric in list : "<<pname<<" -> "<<cmetric<<"\n";
  }
  else {
    //initial bases calc 
    //Compute metric for process arranged like q...qb...g
    msg_Debugging() << name << "\n";
    cmetric=CMetric_Base::GetCM(CMetric_Key(name,tmp));
    if (cmetric==nullptr) THROW(not_implemented,"No metric for "+name);
    msg_Debugging()<<"Metric for '"<<pname<<"' is "<<cmetric<<"\n";
    m_cmetrics.insert(make_pair(pname,cmetric));
  }

  m_ordered_ids_np1 = vector<size_t>(tmp->Legs().size());
  for(size_t i=0; i<tmp->Legs().size(); i++) {
    m_ordered_ids_np1.at(i) = tmp->Leg(cmetric->Map(i))->Id();
  }


  
  p_cmetric_np1=cmetric;

  //TODO: choose momenta at random, this probably needs much improvement still
  size_t leg_emit, leg_soft, leg_spect;

  for(size_t i=p_ampl_np1->NIn();i<p_ampl_np1->Legs().size();i++) {
    if(p_ampl_np1->Leg(i)->Flav().Kfcode() == 21) {
      leg_soft = i;
      if(i==p_ampl_np1->NIn()) {
        leg_emit = i+1;
        leg_spect = i+2;
      }
      else if(i==p_ampl_np1->Legs().size()-1) {
        leg_emit = i-2;
        leg_spect = i-1;
      }
      else {
        leg_emit = i-1;
        leg_spect = i+1;
      }
      break;
    }
  }
  
  kf_code new_fl = new_flavour(p_ampl_np1->Leg(leg_emit)->Flav().Kfcode(),
                               p_ampl_np1->Leg(leg_soft)->Flav().Kfcode());

  p_emit = p_ampl_np1->Leg(leg_emit);
  p_soft = p_ampl_np1->Leg(leg_soft);
  p_spect = p_ampl_np1->Leg(leg_spect);

  
  Vec4D p1 = p_emit->Mom();
  Vec4D p2 = p_soft->Mom();
  Vec4D p3 = p_spect->Mom();
  
  double y = p1*p2/(p1*p2+p1*p3+p2*p3);

  
  p_ampl_n = p_ampl_np1->Copy();

  p_ampl_n->SetOrderQCD(p_ampl_np1->OrderQCD()-1);
  p_ampl_n->CombineLegs(p_ampl_n->Leg(leg_emit),p_ampl_n->Leg(leg_soft),
                      Flavour(new_fl));
  // p_ampl_n->Leg(leg1)->SetMom(p_ampl_n->Leg(leg1)->Mom()-y/(1-y)*p_ampl_np1->Leg(leg3)->Mom());
  // p_ampl_n->Leg(leg3)->SetMom(1./(1.-y) * p_ampl_n->Leg(leg3)->Mom());
  for(size_t i=0;i<p_ampl_n->Legs().size();i++) {
    msg_Debugging()<<"Set leg "<<i<<" "<<"\n";
    if(i==leg_soft) {
      p_ampl_n->Leg(i)->SetMom(p_emit->Mom()+p_soft->Mom()-y/(1-y)*p_spect->Mom());
      p_emit_n = p_ampl_n->Leg(i);
      p_ampl_n->Leg(i+1)->SetMom(1./(1.-y) * p_spect->Mom());
      p_spect_n = p_ampl_n->Leg(i+1);
    }
    p_ampl_n->Leg(i)->SetId(pow(2,i));
    msg_Debugging()<<*p_ampl_n<<"\n";
  }
  msg_Debugging()<<"New amplitude for n-parton process: "<<*p_ampl_n<<"\n";

  tmp = p_ampl_n->Copy();
  tmp->SetNIn(0);
  Process_Base::SortFlavours(tmp);


  
  string pname_n = Process_Base::GenerateName(tmp);

  msg_Debugging()<<"new process: "<<pname_n<<"\n";
  msg_Debugging()<<*tmp<<"\n";

  n_g=0;
  n_q=0;
  n_aq=0;
   
  for (size_t i=0;i<tmp->Legs().size();i++) {
    Flavour flav = tmp->Leg(i)->Flav();
    if (i<2) flav=flav.Bar();
    if (flav==Flavour(kf_gluon)) n_g++;
    if (flav.IsQuark() && !flav.IsAnti()) n_q++;
    if (flav.IsQuark() && flav.IsAnti()) n_aq++;
  }


  CMetric_Base* cmetric_n;
  CMetric_Map::const_iterator CMiter_n=m_cmetrics.find(pname_n);
  string name_n= ToString(n_g)+"_G_"+ToString(n_q)+"_Q_"+ToString(n_aq)+"_AQ";
  msg_Debugging()<<"Found process "<<name_n<<"\n";
  if (CMiter_n!=m_cmetrics.end()) {
    cmetric_n=CMiter_n->second;
    msg_Debugging()<<" found metric in list : "<<pname_n<<" -> "<<cmetric_n<<"\n";
  }
  else {
    //initial bases calc 
    //Compute metric for process arranged like q...qb...g
    cmetric_n=CMetric_Base::GetCM(CMetric_Key(name_n,tmp));
    if (cmetric_n==nullptr) THROW(not_implemented,"No metric for "+name_n);
    msg_Debugging()<<"Metric for '"<<pname_n<<"' is "<<cmetric_n<<"\n";
    m_cmetrics.insert(make_pair(pname_n,cmetric_n));
  }

  m_ordered_ids = vector<size_t>(tmp->Legs().size());
  for(size_t i=0; i<tmp->Legs().size(); i++) {
    m_ordered_ids.at(i) = tmp->Leg(cmetric_n->Map(i))->Id();
  }

  p_cmetric_n = cmetric_n;
  tmp->Delete();
  return true;
} 
  

double RRatios::CplFac(const ATOOLS::Flavour &/*fli*/,const ATOOLS::Flavour &/*flj*/,
                       const ATOOLS::Flavour &/*flk*/,const int /*type*/,
                       const int /*cpl*/,const double &/*mu2*/) const
{
  THROW(not_implemented,"");
  return -1.0;
}




DECLARE_GETTER(RRatios,"RRatios",Shower_Base,Shower_Key);

Shower_Base *Getter<Shower_Base,Shower_Key,RRatios>::
operator()(const Shower_Key &key) const
{
  return new RRatios(key.p_isr,key.p_model);
}

void Getter<Shower_Base,Shower_Key,RRatios>::
PrintInfo(std::ostream &str,const size_t /*width*/) const
{ 
  str<<"The Colorful Resummation"; 
}
