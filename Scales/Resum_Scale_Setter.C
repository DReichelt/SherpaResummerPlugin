#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

#include "Analysis/Observable_Base.H"

using namespace PHASIC;
using namespace ATOOLS;

namespace PHASIC {

  class ResumObs_Scale_Setter: public Scale_Setter_Base {
  protected:

  public:

    ResumObs_Scale_Setter(const Scale_Setter_Arguments &args) :
      Scale_Setter_Base(args)
    {
      m_scale.resize(3); // by default three scales: fac, ren, res
                         // but you can add more if you need for COUPLINGS
      SetCouplings(); // the default value of COUPLINGS is "Alpha_QCD 1", i.e.
                      // m_scale[1] is used for running alpha_s
                      // (counting starts at zero!)
      const std::string& scale = args.m_scale;
      const std::string& name = "ResumObs{";
      size_t start = scale.find(name);
      if(start == std::string::npos) THROW(fatal_error,"Invalid scale '"+args.m_scale+"'");
      start += name.size();
      size_t end = scale.find("}");
      if(end == std::string::npos) THROW(fatal_error,"Invalid scale '"+args.m_scale+"'");
      std::string Obs = scale.substr(start,end-start);
      p_Obs = RESUM::Observable_Base::Ptr(RESUM::Observable_Getter::GetObject(Obs, RESUM::Observable_Key(Obs)));
      if(p_Obs == nullptr) THROW(fatal_error,"Observable not found '"+Obs+"'");
    }

    double Calculate(const std::vector<ATOOLS::Vec4D> &p,
		     const size_t &mode)
    {
      std::vector<ATOOLS::Flavour> dummy(p.size(),Flavour(21));
      // double a = p_Obs->Parameters(p, dummy, 0).m_a;
      double Q2 = (p[0]+p[1]).Abs2(); 
      double muF= Q2;
      double muR= p_Obs->Value(p,dummy)*Q2;
      double muQ= p_Obs->Value(p,dummy)*Q2;

      m_scale[stp::fac] = muF;
      m_scale[stp::ren] = muR;
      m_scale[stp::res] = muQ;

      // Switch on debugging output for this class with:
      // Sherpa "OUTPUT=2[ResumObs_Scale_Setter|15]"
      DEBUG_FUNC("Calculated scales:");
      DEBUG_VAR(m_scale[stp::fac]);
      DEBUG_VAR(m_scale[stp::ren]);
      DEBUG_VAR(m_scale[stp::res]);

      return m_scale[stp::fac];
    }
 private:
  RESUM::Observable_Base::Ptr p_Obs;
  };

  
}

// Some plugin magic to make it available for SCALES=CUSTOM
DECLARE_GETTER(ResumObs_Scale_Setter,"ResumObs",
	       Scale_Setter_Base,Scale_Setter_Arguments);

Scale_Setter_Base *ATOOLS::Getter
<Scale_Setter_Base,Scale_Setter_Arguments,ResumObs_Scale_Setter>::
operator()(const Scale_Setter_Arguments &args) const
{
  return new ResumObs_Scale_Setter(args);
}

void ATOOLS::Getter<Scale_Setter_Base,Scale_Setter_Arguments,
		    ResumObs_Scale_Setter>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"ResumObs scale scheme";
}
