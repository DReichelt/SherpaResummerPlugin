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
  class JetMass_Base: public Observable_Base {
  public:

    JetMass_Base(const Observable_Key &args): 
      Observable_Base(args), m_algkey(args) {
      m_algkey.m_name = "FJmaxPTjet";
      m_algtag = m_algkey.Name();
      m_algtag += ":"+args.KwArg("R","0.8");
      m_algtag += ":"+args.KwArg("minPT","0");
      m_R = to_type<double>(args.KwArg("R","0.8"));

      // soft drop parameters
      m_zcut = to_type<double>(args.KwArg("zcut","0.0"));
      m_beta = to_type<double>(args.KwArg("beta","0"));
      m_R0 = m_R;
      if(GROOM==0 or m_zcut==0.) m_gmode = GROOM_MODE::NONE;
      else m_gmode = GROOM_MODE::SD;
      DEBUG_FUNC(Name()+" -> "+Tag());
      DEBUG_VAR(m_zcut);
      DEBUG_VAR(m_beta);
      DEBUG_VAR(m_R);
      DEBUG_VAR(m_gmode);
    }

    Obs_Params Parameters
      (const std::vector<ATOOLS::Vec4D>& p,
       const std::vector<ATOOLS::Flavour>& fl,
       const size_t& l) {
      // @TODO: returning a=0 leads to problems, because a_0 is magic
      if(l<2) return {1,0,0,0,0,0};
      // @TODO: we could check here if l corresponds to the largest pT jet,
      // possibly enabling this also for the leading di-/multijet etc.
      // at the moment, assume we are in pp -> Zj, so only one colour charged final
      // state that corresponds to the leading jet
      if(!fl[l].Strong()) return {1,0,0,0,0,0};
      const double a = 1;
      const double b = 1;
      const Poincare cms= {p[0]+p[1]};
      Vec4D pl = p[l];
      Poincare(p[0]+p[1]).Boost(pl);
      const double d = 2.*cosh(pl.Eta())/m_R * (p[0]+p[1]).Abs()/(p[l].PPerp()*m_R);
      const double etamin = log(2.*cosh(pl.Eta())/m_R);
      return Obs_Params(a,b,log(d),0.0,etamin,1);
    }

    std::function<double(double,double&)> FFunction(const std::vector<ATOOLS::Vec4D>& p,
                                                    const std::vector<ATOOLS::Flavour>& fl, 
                                                    const RESUM::Params& params) {
      return FFUNCTION::Additive;
    }

    virtual std::set<size_t> ResumMult() {return {3};}
    virtual size_t ResumQCDorderBorn() {return 1;};
    virtual size_t ResumQCDorderLO() {return 2;}
    virtual size_t ResumQCDorderNLO() {return 3;}


    GROOM_MODE GroomMode(double v, ATOOLS::Cluster_Amplitude* ampl,
                         const size_t &l) {
      if(l<2 or not ampl->Leg(l)->Flav().Strong()) {
        return m_gmode;
      }
      else {
        // return GROOM_MODE::SD_COLL;
        const double tp = GroomTransitionPoint(ampl, l);
        if(v > tp) {
          msg_Debugging()<<"Value v = "<<v<<" above transition point "<<tp<<"\n";
          return GROOM_MODE::NONE;
        }
        else {
          msg_Debugging()<<"Value v = "<<v<<" below transition point "<<tp<<"\n";
          return GROOM_MODE::SD_COLL;        
        }
      }
    }
    
//     double GroomTransitionPoint(ATOOLS::Cluster_Amplitude* ampl,
//                                 const size_t &l) {
//       
//       const double sinth = sqrt(2.0*ampl->Leg(2)->Mom().PPerp2() /
//                                 (ampl->Leg(0)->Mom()*ampl->Leg(1)->Mom()));
//       
//       const Obs_Params para = Observable_Base::Parameters(ampl,l);
//       const double logfac = LogFac(ampl);
//       
//       return pow(2,m_beta)*m_zcut/pow(m_R*sinth,m_beta)*exp(para.m_logdbar)/logfac;
//     }


    double SoftGlobal(ATOOLS::Cluster_Amplitude* ampl, size_t i, size_t j, double scale) override {
      if(ampl->Leg(i)->Flav().Strong() and ampl->Leg(j)->Flav().Strong()) {
        if(i<ampl->NIn() and j<ampl->NIn()) {
          return sqr(m_R)/4.;
        }
        else {
          // return sqr(m_R)*(1./4.)/4.;
          return sqr(m_R)*(1./4.+sqr(m_R)/288.)/4.;
        }
      }
      else {
        return 0;
      }
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
      for(const Vec4D& p: ip) {
        if(p.Nan()) {
          // _rejected++;
          // msg_Error()<<METHOD<<": Nan in momenta, rejecting (returning 0)\n";
          // for(const Vec4D& p: ip) msg_Error()<<p<<"\n";
          // msg_Error()<<"Rejected "<<_rejected<<" events so far.\n";
          return 0;
        }
      }


      if(ip.size() <= nin) return 0;
      msg_Debugging()<<"Start jet mass.\n";
      
      auto alg = algorithms.find(m_algtag);
      if(alg==algorithms.end()) {
        alg = algorithms.insert({m_algtag,GetAlgorithm<double>(m_algkey, ip, fl, nin)}).first;
        msg_Debugging()<<"Found jets.\n";
      }
      else msg_Debugging()<<"Reusing jets found earlier.\n";

      const ATOOLS::Vec4D jet =  alg->second->jetVectors(GROOM);
      const double jmass = jet.Abs2() / sqr(m_R*jet.PPerp());
      return jmass;
    }           

  private:
    double m_R;
    std::string m_pstring;
    Observable_Key m_algkey = {"AlgKey"};
    std::string m_algtag;
  };// end of class Y1_II

}// end of namespace RESUM

using namespace RESUM;

typedef JetMass_Base<0> JetMass;
DECLARE_GETTER(JetMass,"JetMass",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,JetMass>::
operator()(const Parameter_Type &args) const 
{ return new JetMass(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,JetMass>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"JetMass"; }

typedef JetMass_Base<1> SD_JetMass;
DECLARE_GETTER(SD_JetMass,"SD_JetMass",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,SD_JetMass>::
operator()(const Parameter_Type &args) const 
{ return new SD_JetMass(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,SD_JetMass>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"SD_JetMass"; }

