#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "Observables/YN_Durham.H"
#include "Tools/StringTools.H"

#include "Analysis/Observable_Base.H"

using namespace PHASIC;
using namespace ATOOLS;

namespace PHASIC {

  class ResumObs_Scale_Setter_Durham: public Scale_Setter_Base {
  protected:

  public:

    ResumObs_Scale_Setter_Durham(const Scale_Setter_Arguments &args) :
      Scale_Setter_Base(args)
    {
      m_scale.resize(3); // by default three scales: fac, ren, res
                         // but you can add more if you need for COUPLINGS
      SetCouplings(); // the default value of COUPLINGS is "Alpha_QCD 1", i.e.
                      // m_scale[1] is used for running alpha_s
                      // (counting starts at zero!)
      const std::string& scale = args.m_scale;
      const std::string& name = "ResumDurham{";
      size_t start = scale.find(name);
      if(start == std::string::npos) THROW(fatal_error,"Invalid scale '"+args.m_scale+"'");
      start += name.size();
      size_t end = scale.find("}");
      if(end == std::string::npos) THROW(fatal_error,"Invalid scale '"+args.m_scale+"'");
      std::vector<std::string> ags = RESUM::split(scale.substr(start,end-start),",");
      m_n = std::stoi(ags[0]);
      m_murFac = ags.size() > 1 ? std::stod(ags[1]) : 1.;
      m_muqFac = ags.size() > 2 ? std::stod(ags[2]) : 1.;
      m_mufFac = ags.size() > 3 ? std::stod(ags[3]) : 1.;
      std::string Obs = "Y3_Durham";
      p_y = (RESUM::YN_Durham<3>*) (RESUM::Observable_Getter::GetObject(Obs, RESUM::Observable_Key(Obs)));
      p_y4 = (RESUM::YN_Durham<4>*) (RESUM::Observable_Getter::GetObject("Y4_Durham", RESUM::Observable_Key("Y4_Durham")));
      if(p_y == nullptr) THROW(fatal_error,"Observable not found '"+Obs+"'");
    }

    double Calculate(const std::vector<ATOOLS::Vec4D> &p,
		     const size_t &mode)
    {
      std::vector<ATOOLS::Flavour> dummy(p.size(),Flavour(21));
      // double a = p_Obs->Parameters(p, dummy, 0).m_a;
      double Q2 = (p[0]+p[1]).Abs2(); 
      double muF= Q2;
      
      MODEL::Running_AlphaS* as=(MODEL::Running_AlphaS*)MODEL::s_model->GetScalarFunction("alpha_S");
      const double beta0 = as->Beta0(Q2)/(2.*M_PI);
      const double aS = (*as)(Q2);
      const std::vector<double>& scales = p_y->AllValues(p,dummy,2,m_n);
      double prod = 1.;
      for(double y: scales) {
        double lambda = -aS*beta0*log(y);
        if(lambda > 1) msg_Error()<<"Scale below lambda QCD, lambda = "<<lambda<<"\n";
        // msg_Out()<<y<<"\n";
        prod *= (1.-lambda);
      }
      prod = pow(prod,1./(scales.size()+1));
      prod -= 1.;
      prod /= aS*beta0;
      // msg_Out()<<"y3 "<<p_y->Value(p,dummy,2)<<"\n";
      // msg_Out()<<"y4 "<<p_y4->Value(p,dummy,2)<<"\n";
      // msg_Out()<<exp(prod)<<"\n\n";
      double muR= Q2*exp(prod);
      double muQ= scales[0]*Q2;
      // msg_Out()<<sqrt(muQ)<<" "<<sqrt(muR)<<"\n";
      m_scale[stp::fac] = m_mufFac*muF;
      m_scale[stp::ren] = m_murFac*muR;
      m_scale[stp::res] = m_muqFac*muQ;

      // Switch on debugging output for this class with:
      // Sherpa "OUTPUT=2[ResumObs_Scale_Setter|15]"
      DEBUG_FUNC("Calculated scales:");
      DEBUG_VAR(m_scale[stp::fac]);
      DEBUG_VAR(m_scale[stp::ren]);
      DEBUG_VAR(m_scale[stp::res]);

      return m_scale[stp::fac];
    }
  private:
    int m_n = -1;
    double m_murFac = 1.;
    double m_muqFac = 1.;
    double m_mufFac = 1.;
    RESUM::YN_Durham<3>* p_y;
    RESUM::YN_Durham<4>* p_y4;
    RESUM::YN_Durham<5>* p_y5;
  };

  
}

// Some plugin magic to make it available for SCALES=CUSTOM
DECLARE_GETTER(ResumObs_Scale_Setter_Durham,"ResumDurham",
	       Scale_Setter_Base,Scale_Setter_Arguments);

Scale_Setter_Base *ATOOLS::Getter
<Scale_Setter_Base,Scale_Setter_Arguments,ResumObs_Scale_Setter_Durham>::
operator()(const Scale_Setter_Arguments &args) const
{
  return new ResumObs_Scale_Setter_Durham(args);
}

void ATOOLS::Getter<Scale_Setter_Base,Scale_Setter_Arguments,
		    ResumObs_Scale_Setter_Durham>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"ResumObs scale scheme with durham scale.";
}
