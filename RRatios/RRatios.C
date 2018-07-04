#include "RRatios/RRatios.H"

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



RRatios::RRatios(ISR_Handler *const isr,
	     Model_Base *const model):
  Shower_Base("RRatios"), p_ampl(NULL)
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

RRatios::~RRatios()
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

inline kf_code new_flavour(kf_code fl1, kf_code fl2) {
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
  // TODO: this should of course at some point become a loop
  Vec4D emit = p_emit_n->Mom();
  Vec4D spect = p_spect_n->Mom();

  double lambda = 0.001;
  Vec4D soft = lambda*p_soft->Mom();
  double eps = emit*soft/(spect*(emit-soft));
  Vec4D p1 = emit-soft + eps * spect;
  Vec4D p3 = (1.-eps)*spect;


  msg_Debugging()<<"New momenta:\n"<<p1<<"\n"<<soft<<"\n"<<p3<<"\n";
  
  p_emit->SetMom(p1);
  p_soft->SetMom(soft);
  p_spect->SetMom(p3);

  msg_Debugging()<<*p_ampl_np1<<"\n";
  
  MatrixD metric_np1 = p_cmetric_np1->CMetric();
  MatrixD metric_n = p_cmetric_n->CMetric();

  metric_np1.print(msg_Out());
  metric_n.print(msg_Out());
  
  unsigned int dim_np1 = metric_np1.dim();  
  unsigned int dim_n = metric_n.dim();

  
  MatrixD H_np1 = MatrixC(p_comix->ComputeHardMatrix(p_ampl_np1,p_cmetric_np1->Perms()),
                          dim_np1, 0).real();
  MatrixD H_n = MatrixC(p_comix->ComputeHardMatrix(p_ampl_n,p_cmetric_n->Perms()),
                        dim_n, 0).real();
 
  // TODO: this currently only exists to create debugging output
  Vec4D_Vector moms(p_ampl->Legs().size());
  Flavour_Vector flavs(p_ampl->Legs().size());
  for (size_t i(0);i<p_ampl->Legs().size();++i) {
    moms[i]=i<p_ampl->NIn()?-p_ampl->Leg(i)->Mom():p_ampl->Leg(i)->Mom();
    flavs[i]=i<p_ampl->NIn()?p_ampl->Leg(i)->Flav().Bar():p_ampl->Leg(i)->Flav();
  }

  
  double s = 2.*moms[0]*moms[1];
  double t = -2.*moms[0]*moms[2];
  double u = -2.*moms[0]*moms[3];
  //

  
  MatrixD cH_np1 = metric_np1*H_np1;
  MatrixD cH_n = metric_n*H_n;
  
  std::cout << "c*H:" << std::endl;
  cH_np1.print(std::cout);
  
  

  std::cout << "c_n*H_n:" << std::endl;
  cH_n.print(std::cout);
  
  
  double TrcH = cH_np1.trace();
  double TrcH_n = cH_n.trace();
  
  std::cout<< "Tr c*H = " << TrcH/1536/0.3302891295379082/0.3302891295379082  * 32 << "\n";
  std::cout<< "Tr c_n*H_n = " << TrcH_n/512/0.3611575592573076/0.3611575592573076 * 16 << "\n";


  
  std::vector<MatrixD> Tprods(p_cmetric_n->Tprods().size());
  for(size_t i=0; i<p_cmetric_n->Tprods().size(); i++) Tprods.at(i) = p_cmetric_n->Tprods().at(i);



  

  MatrixD Gamma(Tprods.at(0).numElements(), Tprods.at(0).dim());
  MatrixD SumTs(Tprods.at(0).numElements(), Tprods.at(0).dim());
  MatrixD Imetric = p_cmetric_n->Imetric();

  // MatrixD Imetric = {{{268823., 268822., 268822., -194174., -194174., -194174.},
  //                     {268822., 268823., 268822., -194174., -194174., -194174.},
  //                     {268822., 268822., 268823., -194174., -194174., -194174.},
  //                     {-194174., -194174., -194174., 140256., 140254., 140254.},
  //                     {-194174., -194174., -194174., 140254., 140256., 140254.},
  //                     {-194174., -194174., -194174., 140254., 140254., 140256.}}};

  int r = 0;
  int f = 0;

  msg_Out()<<*p_ampl_np1<<"\n";
  msg_Out()<<soft<<"\n";
  msg_Out()<<*p_ampl_n<<"\n";
  
  msg_Out()<<"\n";
  metric_n.print(msg_Out());
  msg_Out()<<"\n";
  Imetric.print(msg_Out());
  msg_Out()<<"\n";
  (metric_n*Imetric).print(msg_Out());
  msg_Out()<<"\n";
  size_t i=0;
  for(size_t f=0; f<p_ampl_n->Legs().size(); f++){ 
    for(size_t r=f+1; r<p_ampl_n->Legs().size(); r++) {
      msg_Out()<<f<<" "<<r<<"\n";
      Vec4D pt = p_ampl_n->Leg(f)->Mom();
      Vec4D pr = p_ampl_n->Leg(r)->Mom();
      if(f < p_ampl_n->NIn()) pt *= -1;
      if(r < p_ampl_n->NIn()) pr *= -1;
      double eikonal = pt*pr / ((pt*soft)*(pr*soft));
      (Tprods.at(i)).print(std::cout);
      msg_Out()<<pt<<"\n";
      msg_Out()<<pr<<"\n";
      msg_Out()<<eikonal<<"\n";
      Gamma += eikonal*Tprods.at(i);
      SumTs += Tprods.at(i);
      i++;
    }
  }
  msg_Out()<<"Sum Cij\n";
  SumTs.print(msg_Out());
  msg_Out()<<"\n";
  ((2.*3.*MatrixD::diagonal(1,0,SumTs.dim())+SumTs)).print(msg_Out());
  msg_Out()<<"\n";
  MatrixD chGamma = metric_n*H_n*Gamma;
  chGamma.print(msg_Out());
  double TrcHG = 0.118/M_PI * chGamma.trace();
  std::cout<< "Tr c_n*H_n = " << TrcHG/512/0.3611575592573076/0.3611575592573076 * 16 << "\n";
  exit(1);
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

void RRatios::CleanUp()
{
  p_comix->Reset();
  if (p_ampl) p_ampl->Delete();
  p_ampl=NULL;  
}

Cluster_Definitions_Base *RRatios::GetClusterDefinitions()
{
  return NULL;
}


bool RRatios::PrepareShower
(Cluster_Amplitude* ampl,const bool &soft)
{
  DEBUG_FUNC(this);
  DEBUG_VAR(ampl->Proc<Process_Base>());


  p_ampl_np1=ampl->Copy();
  p_ampl=p_ampl_np1;
  //p_ampl = p_ampl_np1;

  // if we want we can reset momenta here
  Vec4D emit = {100,50,50,0};
  Vec4D spect = {100,50,-50,0};

  double lambda = 0.001;
  Vec4D s = {lambda,lambda,0,0};
  double eps = emit*s/(spect*(emit-s));
  emit += -s + eps * s;
  spect *= (1.-eps);
  SetMomenta(p_ampl_np1, {{100,0,0,100},{100,0,0,-100},emit,s,spect});
  msg_Out()<<*p_ampl_np1<<"\n";

  p_ampl_np1->SetNIn(0);  
  std::string pname=Process_Base::GenerateName(p_ampl_np1);

  Process_Base::SortFlavours(p_ampl_np1);
  CMetric_Base* cmetric;

  n_g=0;
  n_q=0;
  n_aq=0;
   
  for (size_t i=0;i<p_ampl_np1->Legs().size();i++) {
    Flavour flav = p_ampl_np1->Leg(i)->Flav();
    if (i<2) flav=flav.Bar();
    if (flav==Flavour(kf_gluon)) n_g++;
    if (flav.IsQuark() && !flav.IsAnti()) n_q++;
    if (flav.IsQuark() && flav.IsAnti()) n_aq++;
  }

  //number of color singlets
  color_sings = p_ampl_np1->Legs().size() - (n_g + n_aq + n_q);
   
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
    cmetric=CMetric_Base::GetCM(CMetric_Key(name,p_ampl_np1));
    if (cmetric==NULL) THROW(not_implemented,"No metric for "+name);
    msg_Debugging()<<"Metric for '"<<pname<<"' is "<<cmetric<<"\n";
    m_cmetrics.insert(make_pair(pname,cmetric));
  }


  
  p_cmetric_np1=cmetric;

  p_cmetric = p_cmetric_np1;

  p_ampl_np1=ampl->Copy();

  //TODO: choose momenta at random
  size_t leg1 = p_ampl_np1->NIn();
  size_t leg2 = leg1+1;
  size_t leg3 = leg2+1;
  kf_code new_fl = new_flavour(p_ampl_np1->Leg(leg1)->Flav().Kfcode(),
                               p_ampl_np1->Leg(leg2)->Flav().Kfcode());

  p_emit = p_ampl_np1->Leg(leg1);
  p_soft = p_ampl_np1->Leg(leg2);
  p_spect = p_ampl_np1->Leg(leg3);
  
  Vec4D p1 = p_emit->Mom();
  Vec4D p2 = p_soft->Mom();
  Vec4D p3 = p_spect->Mom();
  
  double y = p1*p2/(p1*p2+p1*p3+p2*p3);

  
  p_ampl_n = p_ampl_np1->Copy();
  p_ampl_n->SetOrderQCD(p_ampl_np1->OrderQCD()-1);
  p_ampl_n->CombineLegs(p_ampl_n->Leg(leg1),p_ampl_n->Leg(leg2),
                      Flavour(new_fl));
  // p_ampl_n->Leg(leg1)->SetMom(p_ampl_n->Leg(leg1)->Mom()-y/(1-y)*p_ampl_np1->Leg(leg3)->Mom());
  // p_ampl_n->Leg(leg3)->SetMom(1./(1.-y) * p_ampl_n->Leg(leg3)->Mom());
  p_ampl_n->Leg(leg1)->SetMom({100,-50,-50,0});
  p_ampl_n->Leg(leg3)->SetMom({100,50,50,0});
  p_ampl_n->Leg(leg1)->SetId(4);
  p_ampl_n->Leg(leg3)->SetId(8);

  p_emit_n = p_ampl_n->Leg(leg1);
  p_spect_n = p_ampl_n->Leg(leg3);

  ATOOLS::Cluster_Amplitude* tmp = p_ampl_n->Copy();
  tmp->SetNIn(0);
  Process_Base::SortFlavours(tmp);
  std::string pname_n = Process_Base::GenerateName(tmp);

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
  std::string name_n= ToString(n_g)+"_G_"+ToString(n_q)+"_Q_"+ToString(n_aq)+"_AQ";
  msg_Debugging()<<"Found process "<<name_n<<"\n";
  if (CMiter_n!=m_cmetrics.end()) {
    cmetric_n=CMiter_n->second;
    msg_Debugging()<<" found metric in list : "<<pname_n<<" -> "<<cmetric_n<<std::endl;
  }
  else {
    //initial bases calc 
    //Compute metric for process arranged like q...qb...g
    cmetric_n=CMetric_Base::GetCM(CMetric_Key(name_n,tmp));
    if (cmetric_n==NULL) THROW(not_implemented,"No metric for "+name_n);
    msg_Debugging()<<"Metric for '"<<pname_n<<"' is "<<cmetric_n<<"\n";
    m_cmetrics.insert(make_pair(pname_n,cmetric_n));
  }

  p_cmetric_n = cmetric_n;
  tmp->Delete();
  return true;
} 
  

double RRatios::CplFac(const ATOOLS::Flavour &fli,const ATOOLS::Flavour &flj,
		    const ATOOLS::Flavour &flk,const int type,
		    const int cpl,const double &mu2) const
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
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"The Colorful Resummation"; 
}
