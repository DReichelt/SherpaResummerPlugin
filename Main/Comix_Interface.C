#include "Main/Comix_Interface.H"

#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/ME_Generators.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "COMIX/Main/Single_Process.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_File.H"

using namespace RESUM;
using namespace PHASIC;
using namespace ATOOLS;

namespace RESUM {

  std::ostream &operator<<(std::ostream &s,const Hard_Coefficient &c)
  {
    return s<<"{id="<<c.m_id<<",amp="<<c.m_amps<<"}";
  }

}

Comix_Interface::Comix_Interface(): p_h(nullptr)
{
  m_pmap[nlo_type::lo] = new StringProcess_Map();
}

Comix_Interface::~Comix_Interface()
{
  for (size_t i(0);i<m_procs.size();++i) delete m_procs[i];
  delete m_pmap[nlo_type::lo];
  if (p_h) delete p_h;
}


Hard_Matrix *Comix_Interface::ComputeHardMatrix
(ATOOLS::Cluster_Amplitude *const ampl,
 const std::vector<PHASIC::Idx_Vector>& perms) {
  if(not p_h) {
    COMIX::Single_Process *xs=
      GetProcess(ampl)->Get<COMIX::Single_Process>();
    DEBUG_FUNC(xs->Name());
    msg_Debugging()<<*ampl<<"\n";
    Vec4D_Vector p(ampl->Legs().size());
    for (size_t i(0);i<p.size();++i)
      p[i]=i<ampl->NIn()?-ampl->Leg(i)->Mom():ampl->Leg(i)->Mom();
    xs->GetAmplitude()->SetMomenta(p);
    std::vector<double> s(xs->ScaleSetter(1)->Scales().size(),0.0);
    s[stp::fac]=ampl->MuF2();
    s[stp::ren]=ampl->MuR2();
    s[stp::res]=ampl->MuQ2();
    if (s.size()>stp::size+stp::res) s[stp::size+stp::res]=ampl->KT2();
    xs->SetFixedScale(s);
    xs->ScaleSetter(1)->CalculateScale(p);
    std::vector<Hard_Coefficient> hc;
    for (size_t k(0);k<perms.size();++k) {
      const Idx_Vector &perm(perms[k]);
      msg_Debugging()<<"Permutation "<<perm<<std::endl;
      PHASIC::Int_Vector ci(perm.size(),0), cj(perm.size(),0);
      int fc(ampl->Leg(perm.front())->Flav().StrongCharge());
    
      for (size_t i(0);i<perm.size();++i) {
        size_t cur(perm[i]), next(i+1<perm.size()?perm[i+1]:perm[0]);
        msg_Debugging()<<fc<<" "<<ampl->Leg(cur)->Flav().StrongCharge()<<"\n";
        if (ampl->Leg(cur)->Flav().StrongCharge()==0) {
          continue;
        }
        if (ampl->Leg(cur)->Flav().StrongCharge()==-fc) {
          fc=ampl->Leg(next)->Flav().StrongCharge();
          continue;
        }
        if(ampl->Leg(next)->Flav().Strong()) {
          cj[fc>0?next:cur]=ci[fc>0?cur:next]=Flow::Counter();
        }
        else {
          if(fc!=8) THROW(fatal_error,"Inconsistent colors.");
          cj[perm[0]]=ci[cur]=Flow::Counter();
        }
      }
      double me(xs->GetAmplitude()->Differential(ci,cj,-1));
      std::vector<Spin_Amplitudes> amps;
      std::vector<std::vector<Complex> > cols;
      xs->FillAmplitudes(amps,cols);
      hc.push_back(Hard_Coefficient(perm,amps.front()));
      msg_Debugging()<<"Add "<<hc.back().m_id
                     <<" ( me2 = "<<me<<" ) {\n";
      {
        msg_Indent();
        msg_Debugging()<<"ci = "<<ci<<"\n";
        msg_Debugging()<<"cj = "<<cj<<"\n";
        msg_Debugging()<<hc.back().m_amps;
      }
      msg_Debugging()<<"}\n";
    }
    p_h = new Hard_Matrix(hc.size());
    //if(perms.size()==1){
    //  (*p_h)[0][0]=1; 
    //  return p_h;
    //}
    for (size_t i(0);i<hc.size();++i) {
      p_h->m_id[i]=hc[i].m_id;
      for (size_t j(0);j<hc.size();++j) {
        for (size_t k(0);k<hc[i].m_amps.size();++k) {
          (*p_h)[i][j]+=hc[i].m_amps[k]*std::conj(hc[j].m_amps[k]);
        }
      }
    }
    xs->SetFixedScale(std::vector<double>());
  }
  msg_Debugging()<<*p_h<<"\n";
  return p_h;
}


Process_Base *Comix_Interface::GetProcess
(ATOOLS::Cluster_Amplitude *const ampl)
{
  NLOTypeStringProcessMap_Map *mp=
    ampl->Procs<NLOTypeStringProcessMap_Map>();

  if (not mp) THROW(fatal_error,"Missing NLO process map");
  StringProcess_Map *pm((*mp)[nlo_type::lo]);
  if(not pm) pm = (*mp)[nlo_type::born];
  if (not pm) THROW(fatal_error,"Missing process map");
  Process_Base::SortFlavours(ampl);

  std::string name(Process_Base::GenerateName(ampl));

  StringProcess_Map::const_iterator pit(pm->find(name));
  Process_Base *xs=nullptr;

  if (pit!=pm->end() && pit->second->
      Integrator()->ColorIntegrator()!=nullptr) {
    xs=pit->second;
  }
  if (xs==nullptr) {

    pm=m_pmap[nlo_type::lo];
    if ((pit=pm->find(name))!=pm->end()) xs=pit->second;
    else {
      MakeDir(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process",true);
      My_In_File::OpenDB
	(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Comix/");
      My_In_File::OpenDB
	(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Sherpa/");
      Process_Info pi;
      pi.m_megenerator="Comix";
      for (size_t i(0);i<ampl->NIn();++i) {
	Flavour fl(ampl->Leg(i)->Flav().Bar());
	if (Flavour(kf_jet).Includes(fl)) fl=Flavour(kf_jet);
	pi.m_ii.m_ps.push_back(Subprocess_Info(fl,"",""));
      }
      for (size_t i(ampl->NIn());i<ampl->Legs().size();++i) {
	Flavour fl(ampl->Leg(i)->Flav());
	if (Flavour(kf_jet).Includes(fl)) fl=Flavour(kf_jet);
	pi.m_fi.m_ps.push_back(Subprocess_Info(fl,"",""));
      }
      PHASIC::Process_Base *proc=
	ampl->Proc<Process_Base>()->
	Generator()->Generators()->InitializeProcess(pi,false);
      if (proc==nullptr) {
	My_In_File::CloseDB
	  (rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Comix/");
        My_In_File::CloseDB
	  (rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Sherpa/");
	(*pm)[name]=nullptr;
	return nullptr;
      }
      m_procs.push_back(proc);
      Selector_Key skey(nullptr,nullptr,true);
      proc->SetSelector(skey);
      proc->SetScale
        (Scale_Setter_Arguments
         (MODEL::s_model,rpa->gen.Variable("SCALES"),"Alpha_QCD 1"));
      proc->SetKFactor(KFactor_Setter_Arguments("NO"));
      proc->Get<COMIX::Process_Base>()->Tests();
      proc->FillProcessMap(&m_pmap);
      My_In_File::CloseDB
	(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Comix/");
      if ((pit=pm->find(name))==pm->end()) THROW(fatal_error,"Internal error");
      xs=pit->second;
    }
    if (xs==nullptr) return nullptr;
  }
  return xs;
}

void Comix_Interface::Reset()
{
  if (p_h) delete p_h;
  p_h=nullptr;
}
