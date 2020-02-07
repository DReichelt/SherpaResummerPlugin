#include "Analysis/Observable_Base.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_Limits.H"
#include "FFunction/FFunctions.H"
#include "Tools/StringTools.H"
#include <vector>       
#include <algorithm>       

#include "Observables/Algorithms/GetAlgorithm.H"

using namespace ATOOLS;

namespace RESUM {

  template <int GROOM>
  class JetAngularities_Base: public Observable_Base {
  public:

    JetAngularities_Base(const Observable_Key &args): 
    Observable_Base(args) {
      m_algkey = args;
      m_algkey.m_name = "FJmaxPTjet";
      m_algtag = m_algkey.Name();
      m_algtag += ":"+args.KwArg("R","0.5");
      m_algtag += ":"+args.KwArg("minPT","0");
      m_alpha = to_type<double>(args.KwArg("alpha","2"));
      m_R = to_type<double>(args.KwArg("R","0.5"));
    }

    Obs_Params Parameters
      (const std::vector<ATOOLS::Vec4D>& p,
       const std::vector<ATOOLS::Flavour>& fl,
       const size_t& l) {
      // TODO: this is a dummy
      return Obs_Params(2.0,0.0,0.0,0.0);
    }

    std::function<double(double,double&)> FFunction(const std::vector<ATOOLS::Vec4D>& p,
                                                    const std::vector<ATOOLS::Flavour>& fl, 
                                                    const RESUM::Params& params) {
      // TODO: this is a dummy
      return FFUNCTION::Additive;
    }


    double Value(const std::vector<Vec4D>& moms,
            const std::vector<ATOOLS::Flavour>& flavs,
            const size_t& nin) {
      std::map<std::string, typename Algorithm<double>::Ptr> dummy;
      return Value(moms, flavs, dummy, nin);
    }


    double Value(const std::vector<Vec4D>& ip,
            const std::vector<ATOOLS::Flavour>& fl,
            std::map<std::string, typename Algorithm<double>::Ptr>& algorithms,
            const size_t& nin) {

      if(ip.size() <= nin) return 0;
      msg_Debugging()<<"Start jet angularity.\n";      
      
      auto alg = algorithms.find(m_algtag);
      if(alg==algorithms.end()) {
        alg = algorithms.insert({m_algtag,GetAlgorithm<double>(m_algkey, ip, fl, nin)}).first;
        msg_Debugging()<<"Found jets.\n";
      }
      else msg_Debugging()<<"Reusing jets found earlier.\n";

      std::vector<ATOOLS::Vec4D> constits =  alg->second->apply(ip,GROOM);
      ATOOLS::Vec4D axis(0.,alg->second->jetAxes(GROOM));
      
      
      double lambda = 0.;
      msg_Debugging()<<"Final States:\n";
      for(int i=0; i<ip.size(); i++) msg_Debugging()<<fl[i]<<" "<<ip[i]<<" "<<ip[i].PPerp()<<"\n";
      msg_Debugging()<<"... in jet:\n";
      for (ATOOLS::Vec4D con: constits) {
        // event with grooming , pt is defined without grooming
        double z = con.PPerp() / alg->second->jetVectors(0).PPerp();
        msg_Debugging()<<con<<" -> z = "<<con.PPerp()<< " / "<< alg->second->jetVectors(GROOM).PPerp()<<" = "<<z<<"\n";
        lambda +=  z * pow(con.DR(axis) / m_R, m_alpha);
      }
      msg_Debugging()<<"\n";
      return lambda;
    }           

  private:
    double m_R;
    double m_alpha;
    std::string m_pstring;
    Observable_Key m_algkey = {"AlgKey"};
    std::string m_algtag;
  };// end of class Y1_II

}// end of namespace RESUM

using namespace RESUM;

typedef JetAngularities_Base<0> JetAngularities;
DECLARE_GETTER(JetAngularities,"JetAngularities",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,JetAngularities>::
operator()(const Parameter_Type &args) const 
{ return new JetAngularities(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,JetAngularities>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"JetAngularities"; }

typedef JetAngularities_Base<1> SD_JetAngularities;
DECLARE_GETTER(SD_JetAngularities,"SD_JetAngularities",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,SD_JetAngularities>::
operator()(const Parameter_Type &args) const 
{ return new SD_JetAngularities(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,SD_JetAngularities>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"SD_JetAngularities"; }

