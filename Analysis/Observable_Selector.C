#include "Analysis/Observable_Base.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "PHASIC++/Selectors/Selector.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace ATOOLS;

namespace PHASIC {

  class Observable_Selector: public Selector_Base {
  private:

    RESUM::Observable_Base *p_tau;

    double m_taumin;

  public:

    Observable_Selector(const Selector_Key &key): 
      Selector_Base("Observable_Selector")
    {
      p_proc=key.p_proc;
      m_fl=(Flavour*)&p_proc->Process()->Flavours().front();
      m_sel_log=new Selector_Log(m_name);
      p_tau = RESUM::Observable_Getter::GetObject
	(key[0][0],RESUM::Observable_Key(key[0][0]));
      if (p_tau==NULL) THROW(fatal_error,"Observable not found '"+key[0][0]+"'");
      m_taumin=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][1]));
    }

    ~Observable_Selector() { delete p_tau; }

    bool NoJetTrigger(const ATOOLS::Vec4D_Vector &p)
    {
      return true;
    }

    bool Trigger(const ATOOLS::Vec4D_Vector &p)
    {
      if (p_proc->Process()->Name().find("QCD(S)")) return true;
      double tau=p_tau->Value(&p.front(),m_fl,p.size(),p_proc->NIn());
      bool pass=tau>m_taumin;
      m_sel_log->Hit(1-pass);
      return pass;
    }

    bool JetTrigger(const ATOOLS::Vec4D_Vector &p,
		    ATOOLS::NLO_subevtlist *const sub)
    {
      if (sub->back()->m_i!=sub->back()->m_j) return true;
      double tau=p_tau->Value(&p.front(),sub->back()->p_fl,
			      sub->back()->m_n,p_proc->NIn());
      bool pass=tau>m_taumin;
      m_sel_log->Hit(1-pass);
      return pass;
    }

    void BuildCuts(Cut_Data *) {}

  };

}

using namespace PHASIC;

DECLARE_ND_GETTER(Observable_Selector,"ResumObs",Selector_Base,Selector_Key,true);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,Observable_Selector>::
operator()(const Selector_Key &key) const
{
  if (key.empty() || key[0].size()<2)
    THROW(fatal_error,"Wrong number of parameters");
  return new Observable_Selector(key);
}

void ATOOLS::Getter<Selector_Base,Selector_Key,Observable_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"observable min";
}
