#include "Analysis/Observable_Base.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "PHASIC++/Selectors/Selector.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "Analysis/ChannelAlgorithms/ChannelAlgorithm_Base.H"
#include "Tools/StringTools.H"

using namespace ATOOLS;

namespace PHASIC {

  class Channel_Selector: public Selector_Base {
  private:

    RESUM::ChannelAlgorithm_Base::Ptr p_chAlg = nullptr;

    std::string m_channels;
    enum {SELECT,VETO} m_mode = SELECT;
  public:

    Channel_Selector(const Selector_Key &key): 
      Selector_Base("Channel_Selector")
    {
      p_proc=key.p_proc;
      m_fl=(Flavour*)&p_proc->Process()->Flavours().front();
      m_sel_log=new Selector_Log(m_name);
      RESUM::ChAlg_Key chKey(key[0][0],{key[0].begin(),key[0].end()});
      p_chAlg = RESUM::ChannelAlgorithm_Base::Ptr(RESUM::ChAlg_Getter::GetObject(key[0][0],chKey));
      m_channels = chKey.KwArg("Select","NONE");
      if(m_channels == "NONE") {
        m_mode = VETO;
        m_channels = chKey.KwArg("Veto","");
      }
    }

    ~Channel_Selector() { }

    bool NoJetTrigger(const ATOOLS::Vec4D_Vector &p)
    {
      return true;
    }

    bool Trigger(const ATOOLS::Vec4D_Vector &p)
    {
      const std::string ch = p_chAlg->Channel(p,p_proc->Process()->Flavours(),p_proc->NIn(),true);
      std::map<std::string, typename RESUM::Algorithm<double>::Ptr> algorithms;
      bool pass = m_channels.find(ch) != std::string::npos;
      if(m_mode == VETO) pass = not pass;
      m_sel_log->Hit(1-pass);
      return pass;
    }

    bool JetTrigger(const ATOOLS::Vec4D_Vector &p,
		    ATOOLS::NLO_subevtlist *const sub)
    {
      std::vector<ATOOLS::Flavour> flavs = {sub->back()->p_fl, sub->back()->p_fl+sub->back()->m_n};
      const std::string ch = p_chAlg->Channel(p,flavs,p_proc->NIn(),true);
      bool pass = m_channels.find(ch) != std::string::npos;
      if(m_mode == VETO) pass = not pass;
      m_sel_log->Hit(1-pass);
      return pass;
    }

    void BuildCuts(Cut_Data *) {}

  };

}

using namespace PHASIC;

DECLARE_ND_GETTER(Channel_Selector,"ResumChAlg",Selector_Base,Selector_Key,true);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,Channel_Selector>::
operator()(const Selector_Key &key) const
{
  if (key.empty() || key[0].size()<2)
    THROW(fatal_error,"Wrong number of parameters");
  return new Channel_Selector(key);
}

void ATOOLS::Getter<Selector_Base,Selector_Key,Channel_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Select this channel.";
}
