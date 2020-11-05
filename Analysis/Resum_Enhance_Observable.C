#include "PHASIC++/Enhance/Enhance_Observable_Base.H"
#include "ATOOLS/Math/Algebra_Interpreter.H"
#include "Analysis/Observable_Base.H"

namespace PHASIC {

  class Resum_Enhance_Observable:
    public Enhance_Observable_Base,
    public ATOOLS::Tag_Replacer {
  private:

    ATOOLS::Algebra_Interpreter m_calc;

    const ATOOLS::Vec4D *p_p;

    size_t m_n;

    std::vector<RESUM::Observable_Base *> m_obs;
    std::vector<double> m_obsVals;
    
  public:

    Resum_Enhance_Observable(const Enhance_Arguments &args);

    double operator()(const ATOOLS::Vec4D *p,
		      const ATOOLS::Flavour *fl,
		      const size_t n);

    std::string   ReplaceTags(std::string &expr) const;
    ATOOLS::Term *ReplaceTags(ATOOLS::Term *term) const;

    void AssignId(ATOOLS::Term *term);

  };// end of class Resum_Enhance_Observable

}// end of namespace PHASIC

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "Tools/StringTools.H"

using namespace PHASIC;
using namespace ATOOLS;

DECLARE_GETTER(Resum_Enhance_Observable,"Resum",
	       Enhance_Observable_Base,Enhance_Arguments);

Enhance_Observable_Base *ATOOLS::Getter
<Enhance_Observable_Base,Enhance_Arguments,Resum_Enhance_Observable>::
operator()(const Enhance_Arguments &args) const
{
  return new Resum_Enhance_Observable(args);
}

void ATOOLS::Getter<Enhance_Observable_Base,Enhance_Arguments,Resum_Enhance_Observable>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"variable enhance observable";
}

Resum_Enhance_Observable::Resum_Enhance_Observable
(const Enhance_Arguments &args): Enhance_Observable_Base(args)
{
  std::string arg(args.m_enhance);
  size_t bpos(arg.find("Resum{")), epos(arg.find("}",bpos));
  if (bpos!=0 || epos==std::string::npos)
    THROW(fatal_error,"Invalid input");
  arg = arg.substr(bpos+6,epos-6-bpos);
  m_calc.SetTagReplacer(this);
  bpos = arg.find("Obs[");
  m_obs.clear();
  m_obsVals.clear();
  while(bpos != std::string::npos) {
    epos = arg.find("]",bpos);
    if(epos == std::string::npos) THROW(fatal_error,"Invalid input, missing ']'.");
    const std::string& name = arg.substr(bpos+4,epos-4-bpos);
    std::vector<std::string> ags = RESUM::split(name,",");
    m_obs.push_back(RESUM::Observable_Getter::GetObject(ags[0],RESUM::Observable_Key(ags[0],ags)));
    if(m_obs.back() == nullptr) THROW(fatal_error,"Observable not found '"+name+"'");
    m_calc.AddTag(m_obs.back()->Tag(),"1.0");
    m_obsVals.push_back(1.);
    const std::string& rep = arg.substr(bpos,epos-bpos);
    bpos = arg.find(rep);
    while(bpos != std::string::npos) {
      arg.replace(bpos,rep.size()+1,m_obs.back()->Tag());
      bpos = arg.find(rep);
    }
    bpos = arg.find("Obs[");
  }
  m_calc.Interprete(arg);
}

double Resum_Enhance_Observable::operator()
  (const ATOOLS::Vec4D *p,const ATOOLS::Flavour *fl,const size_t n)
{
  std::vector<ATOOLS::Vec4D> moms = {p,p+n};
  std::vector<ATOOLS::Flavour> flavs = {fl,fl+n};
  for(size_t i=0; i<m_obs.size(); i++) {
    std::map<std::string, typename RESUM::Algorithm<double>::Ptr> algorithms;
    if(m_obs[i]->VetoEvent(moms,flavs,algorithms,p_proc->NIn())) {
      return 0;
    }
    else {
      m_obsVals[i] = m_obs[i]->Value(moms,flavs,algorithms,p_proc->NIn());
    }
  }
  return m_calc.Calculate()->Get<double>();
}

std::string Resum_Enhance_Observable::ReplaceTags(std::string &expr) const
{
  return m_calc.ReplaceTags(expr);
}

Term *Resum_Enhance_Observable::ReplaceTags(Term *term) const
{
  if (term->Id()>=100) {
    if (term->Id()-100>=m_obsVals.size())
      THROW(fatal_error,"Obs index too large");
    term->Set(m_obsVals[term->Id()-100]);
  }
  return term;
}

void Resum_Enhance_Observable::AssignId(Term *term)
{
  for(size_t i=0; i<m_obs.size(); i++) {
    if(term->Tag()==m_obs[i]->Name()) {
      term->SetId(100+i);
      break;
    }
  }
}
