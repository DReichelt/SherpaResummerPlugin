#include "RRatios/RRatios.H"

#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Data_Reader.H"

// temporary yoda dependence
#include "YODA/Scatter2D.h"
#include "YODA/WriterYODA.h" 


using namespace RESUM;
using namespace PDF;
using namespace ATOOLS;
using namespace MODEL;
using PHASIC::Process_Base;
using std::vector;
using std::string;

RRatios::RRatios(ISR_Handler *const /*isr*/,
                 Model_Base *const model):
  Shower_Base("RRatios")
{
  p_as=(Running_AlphaS*)model->GetScalarFunction("alpha_S");
  
  Data_Reader read(" ",";","#","=");
  m_amode=read.GetValue<int>("RESUM_MODE",0);
  rpa->gen.SetVariable("SCALES", read.GetValue<string>("SCALES", "VAR{sqr(91.188)}"));
  if (rpa->gen.Variable("SHOWER_GENERATOR")=="")
    rpa->gen.SetVariable("SHOWER_GENERATOR",ToString(this));
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
  if(skip) return 1;
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
  const size_t dim_np1 = metric_np1.dim();  
  const size_t dim_n = p_cmetric_n->CMetric().size();

  // YODA::Scatter2D::Ptr plot= std::make_shared<YODA::Scatter2D>("/line/line","/line/line");
  YODA::Scatter2D plot("/line/line","/line/line");
  double lambda = 0.99;

  const MatrixD& pref_np1 = p_cmetric_np1->PrefMatrix();
  const MatrixD& pref_n = p_cmetric_n->PrefMatrix(); 
  
  for(double cut=lambda; cut>1e-10; cut*=lambda) {
    const Vec4D& soft = lambda*p_soft->Mom();  
    const double eps = emit*soft/(spect*(emit-soft));
    const Vec4D& p1 = emit-soft + eps * spect;
    const Vec4D& p3 = (1.-eps)*spect;


    msg_Debugging()<<"New momenta:\n"<<p1<<"\n"<<soft<<"\n"<<p3<<"\n\n";
  
    p_emit->SetMom(p1);
    p_soft->SetMom(soft);
    p_spect->SetMom(p3);
  
    MatrixD H_np1 = MatrixC(m_comix.ComputeHardMatrix(p_ampl_np1,
                                                      p_cmetric_np1->Perms()),
                            dim_np1, 0).real();
    H_np1.data() *= pref_np1.data();

    m_comix.Reset();
    MatrixD H_n = MatrixC(m_comix.ComputeHardMatrix(p_ampl_n,p_cmetric_n->Perms()),
                          dim_n, 0).real();
    H_n.data() *= pref_n.data();
    m_comix.Reset(); //TODO: do I need these resets?

    msg_Debugging().precision(20);
    msg_Debugging()<<*p_ampl_np1<<"\n";
    msg_Debugging()<<*p_ampl_n<<"\n";

    
    // msg_Debugging() << "Tr(c_n+1 * H_n+1) = "<< Trace(metric_np1,H_np1)/1536/0.3302891295379082/0.3302891295379082  * 32<<"\n";
    /* msg_Debugging() << "Tr(c_n+1 * H_n+1) = "<< Trace(metric_np1,H_np1)/64./4./0.375/0.375 * 8<<"\n"; */
    // msg_Debugging() << "Tr(c_n * H_n) = "<< Trace(metric_n*H_n)/512/0.3611575592573076/0.3611575592573076 * 16<<"\n\n";
    /* msg_Debugging() << "Tr(c_n * H_n) = "<< Trace(metric_n*H_n)/64./4./0.43301270189221935/0.43301270189221935 * 4<<"\n"; */
    //exit(1);
    // msg_Debugging()<< "Tr c_n+1 * H _n+1  = " << TrcH/1536/0.3302891295379082/0.3302891295379082  * 32 << "\n";
    // // msg_Out()<< "Tr c_n * H_n = " << TrcH_n/512/0.3611575592573076/0.3611575592573076 * 16 << "\n";
    // msg_Debugging()<< "Tr c_n * H_n = " << TrcH_n/4/g/g/g/g<< "\n";

    msg_Debugging() << p_ampl_np1->Leg(3)->Mom().Mass()<<"\n";
    msg_Debugging() << "Tr(c_n+1 * H_n+1) = "<< Trace(metric_np1,H_np1)/64./4.<<"\n";
    msg_Debugging() << "Tr(c_n * H_n) = "<< Trace(metric_n,H_n)/64./4.<<"\n";

    std::vector<MatrixD> Tprods(p_cmetric_n->Tprods().size());
    for(size_t i=0; i<p_cmetric_n->Tprods().size(); i++) Tprods.at(i) = p_cmetric_n->Tprods().at(i);



  

    MatrixD Gamma(Tprods.at(0).dim(), Tprods.at(0).dim());

    msg_Debugging()<<"Amplitude for n+1:"<<*p_ampl_np1<<"\n";
    msg_Debugging()<<"Amplitude for n:"<<*p_ampl_n<<"\n";
  
    size_t i=0;
    for(size_t t=0; t<p_ampl_n->Legs().size(); t++){ 
      for(size_t r=t+1; r<p_ampl_n->Legs().size(); r++) {
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
        msg_Debugging()<<"T-product: \n"<<Tprods.at(i)<<"\n";
        msg_Debugging()<<"p_t = "<<pt<<"\n";
        msg_Debugging()<<"p_r = "<<pr<<"\n";
        msg_Debugging()<<"p_soft = "<<soft<<"\n";
        msg_Debugging()<<"eikonal = "<<eikonal<<"\n\n";
        // TODO: urgent!! why is this a minus and not a plus?????????
        Gamma -= eikonal*Tprods.at(i);
        i++;
      }
    }
    msg_Debugging()<<"Gamma:\n"<<Gamma<<"\n";
    
    double g =  sqrt(4.*M_PI*0.118);
    // double TrcH_np1 = Trace(metric_np1,H_np1) /0.3302891295379082/0.3302891295379082  * 32;

    double TrcH_np1 = Trace(metric_np1,H_np1);///0.375/0.375 * 8;

    
    // we actually want H*metric*colour-change-matrices, but our Tproducts are
    // inverse_metric*colour-change-matrix, so this is correct
    // double TrHG = Trace(H_n,Gamma) /0.3611575592573076/0.3611575592573076 * 16;
    double TrHG = Trace(H_n,Gamma);///0.43301270189221935/0.43301270189221935 * 4;

    
    //msg_Out()<< "Tr c_n*H_n = " << g*g* TrcHG/512/0.3611575592573076/0.3611575592573076 * 16 * 4./6. << "\n";
    // msg_Out()<< "Tr c_n*H_n = " << TrcHG/4./pow(g,4)<< "\n";
    plot.addPoint(cut, (g*g* TrHG) / (TrcH_np1));
    msg_Debugging()<<g*g*TrHG/64./4./2.<<"\n";
    msg_Debugging()<<TrcH_np1/64./4./2.<<"\n";
    msg_Out()<<*p_ampl_np1<<"\n";
    msg_Out()<<*p_ampl_n<<"\n";
    msg_Out()<<(g*g* TrHG)/ (TrcH_np1)<<" "<<cut<<" "<<(g*g* TrHG)<<" "<<(TrcH_np1)<<"\n";
  }
  // YODA::WriterYODA::write(std::to_string(m_count)+".yoda",plot);
  msg_Out()<<"Write out yoda.\n";
  m_count++;
  //exit(1);
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


bool RRatios::PrepareShower
(Cluster_Amplitude* ampl,const bool & /*soft*/)
{
  DEBUG_FUNC(this);
  DEBUG_VAR(ampl->Proc<Process_Base>());

  // if(ampl->Legs().size() != 5) {
  //   skip = true;
  //   return skip;
  // }
  // else {
  //   skip = false;
  // }
  p_ampl_np1=ampl->Copy();


  // if we want we can reset momenta here
  /* Vec4D emit = {7000.0, 2163.1189606246307, -6657.395614066076, 0.}; */
  /* Vec4D spect = {7000.0, -2163.11896062463, 6657.395614066076, 0.}; */

  /* double lambda = 0.001; */
  /* //Vec4D s = {lambda,lambda,0,0}; */
  /* Vec4D s = {140.0, 126.13564150633869, 60.743723476458136, 0.0}; */
  /* double eps = emit*s/(spect*(emit-s)); */
  /* Vec4D emit_rec = emit -s + eps * spect; */
  /* Vec4D spect_rec = spect*(1.-eps); */
  /* SetMomenta(p_ampl_np1, {{7000.,0.,0.,7000.},{7000.,0.,0.,-7000.},s,emit_rec,spect_rec}); */
  /* msg_Out()<<*p_ampl_np1<<"\n"; */

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
