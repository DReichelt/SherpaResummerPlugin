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
      Observable_Base(args), m_algkey(args) {
      m_algkey.m_name = "FJmaxPTjet";
      m_alpha = to_type<double>(args.KwArg("alpha","2"));
      m_R = to_type<double>(args.KwArg("R","0.8"));
      m_WTA = std::set<std::string>({"no","NO","n","N","0"}).count(args.KwArg("WTA","no")) == 0;

      // soft drop parameters
      m_zcut = to_type<double>(args.KwArg("zcut","0.0"));
      m_beta = to_type<double>(args.KwArg("beta","0"));
      m_R0 = m_R;

      m_algtag = m_algkey.Name();
      m_algtag += ":"+args.KwArg("R",std::to_string(m_R));
      m_algtag += ":"+args.KwArg("minPT","0");
      m_algtag += ":"+args.KwArg("zcut",std::to_string(m_zcut));
      m_algtag += ":"+args.KwArg("beta",std::to_string(m_beta));

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
      // updat: do that now, still assume that all final partons correspond to 
      // individual jets
      if(!fl[l].Strong()) return {1,0,0,0,0,0};
      const double pTl = p[l].PPerp();
      for(size_t i=2; i<p.size(); i++) {
        if(fl[i].Strong() and p[i].PPerp() > pTl) return {1,0,0,0,0,0};
      }
      const double a = 1;
      const double b = m_alpha-1;
      const Poincare cms= {p[0]+p[1]};
      Vec4D pl = p[l];
      Poincare(p[0]+p[1]).Boost(pl);
      const double d = pow(2.*cosh(pl.Eta())/m_R,m_alpha-1) * (p[0]+p[1]).Abs()/(p[l].PPerp()*m_R);
      const double etamin = log(2.*cosh(pl.Eta())/m_R);
      return Obs_Params(a,b,log(d),0.0,etamin,1);
    }

    std::function<double(double,double&)> FFunction(const std::vector<ATOOLS::Vec4D>& p,
                                                    const std::vector<ATOOLS::Flavour>& fl, 
                                                    const RESUM::Params& params) {
      return FFUNCTION::Additive;
    }

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
      // See arXiv:1207.1640
      // Note difference of factor 1./4. in definition of f_{l.a.} vs. T(L).
      if(not ampl) THROW(fatal_error, "No amplitude provided.");
      DEBUG_FUNC(*ampl);
      if(ampl->Leg(i)->Flav().Strong() and ampl->Leg(j)->Flav().Strong()) {
        if(i<ampl->NIn() and j<ampl->NIn()) {
          msg_Debugging()<<"Case: Initial + Initial dipole.\n";
          return sqr(m_R)/4.;
        }
        else {
          msg_Debugging()<<"At least one of i,j is final state, so figure out if it is the leading jet.\n";
          int lead = -1;
          double pTlead = -1;
          for(size_t k=ampl->NIn(); k<ampl->Legs().size(); k++) {
            if(ampl->Leg(k)->Flav().Strong() and ampl->Leg(k)->Mom().PPerp() > pTlead) {
              lead = k;
              pTlead = ampl->Leg(k)->Mom().PPerp();
            }
          }
          if(lead<0 or pTlead <0) 
            THROW(fatal_error, "Did not find leading jet.");

          if(i<ampl->NIn() or j<ampl->NIn()) {
            const size_t initial = i<ampl->NIn() ? i : j; 
            const size_t final = i<ampl->NIn() ? j : i;
            if(final==lead) {
              msg_Debugging()<<"Case: Leading + Initial dipole.\n";
              return sqr(m_R)*(1./4.)/4.;
              // return sqr(m_R)*(1./4.+sqr(m_R)/288.)/4.;
            }
            else {
              msg_Debugging()<<"Case: Recoil + Initial dipole.\n";
              // TODO: figure out which order of the signs is correct
              const double dY = (initial==0 ? 1.:-1.) * ampl->Leg(final)->Mom().DY(ampl->Leg(lead)->Mom());
              return 1./8. * exp(dY)/(1.+cosh(dY));
            }
          }
          else {
            // both final state jets
            if(lead==i or lead==j) {
              msg_Debugging()<<"Case: Leading + Recoil dipole.\n";
              return 1./16. * sqr(tanh(ampl->Leg(i)->Mom().DY(ampl->Leg(j)->Mom())));
            }
            else {
              THROW(fatal_error, "Gamma not implemented for multijets.");
            }
          }
        }
      }
      else {
        msg_Debugging()<<"Case: Colour singlet.\n";
        return 0;
      }
    }

    double SoftNonGlobal(ATOOLS::Cluster_Amplitude* ampl, size_t i, size_t j, double scale) override {
      if(ampl->Leg(i)->Flav().Strong() and ampl->Leg(j)->Flav().Strong()) {
        if(i<ampl->NIn() and j<ampl->NIn()) {
          return 4.*pow(m_R,2)*(1.17-log(2.*m_R));
        }
        else {
          return pow(M_PI,2)/3.;
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

#ifdef USING_FJCONTRIB
    double Value(const fastjet::PseudoJet &jet) const{
      // get the jet constituents
      std::vector<fastjet::PseudoJet> constits = jet.constituents();
      if (constits.size() == 0) return 0.0;
    
      // get the reference axis
      fastjet::PseudoJet reference_axis = _get_reference_axis(jet);

      // do the actual coputation
      double numerator = 0.0, denominator = 0.0;
      for (const auto &c : constits){
        double pt = c.pt();
        // Note: better compute (dist^2)^(alpha/2) to avoid an extra square root
        numerator   += pt * pow(c.squared_distance(reference_axis), 0.5*m_alpha);
        denominator += pt;
      }
      return numerator/(denominator*pow(m_R, m_alpha));
    }

    fastjet::PseudoJet _get_reference_axis(const fastjet::PseudoJet &jet) const{
      if (m_alpha>1) return jet;

      fastjet::Recluster recluster(fastjet::JetDefinition(fastjet::antikt_algorithm, fastjet::JetDefinition::max_allowable_R,  fastjet::WTA_pt_scheme));
      return recluster(jet);
    }
#endif

    double Value(const std::vector<Vec4D>& ip,
            const std::vector<ATOOLS::Flavour>& fl,
            std::map<std::string, typename Algorithm<double>::Ptr>& algorithms,
            const size_t& nin) override {

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
      msg_Debugging()<<"Start jet angularity.\n";
      
      auto alg = algorithms.find(m_algtag);
      msg_Debugging()<<"Searching algs...\n";
      if(alg==algorithms.end()) {
        msg_Debugging()<<"Known: \n";
        for(auto alg: algorithms) msg_Debugging()<<alg.first<<"\n";
        msg_Debugging()<<"Did not find "<<m_algtag<<".\n";
        alg = algorithms.insert({m_algtag,GetAlgorithm<double>(m_algkey, ip, fl, nin)}).first;
        msg_Debugging()<<"Found jets.\n";
      }
      else {
        msg_Debugging()<<"Reusing jets found earlier.\n";
      }

#ifdef USING_FJCONTRIB
      const double lambdaFJ = Value(GROOM ? std::dynamic_pointer_cast<FJmaxPTjet>(alg->second)->SDLeadJet() :
                                  std::dynamic_pointer_cast<FJmaxPTjet>(alg->second)->LeadJet());      
      rpa->gen.SetVariable(m_tag, std::to_string(lambdaFJ));
      return lambdaFJ;
#else
      const double lambdaFJ = -1;
#endif

      const std::vector<ATOOLS::Vec4D> constits =  alg->second->apply(ip,GROOM);
      double lambda = 0.;
      if(constits.size() <= 1) {
        msg_Debugging()<<"Jet with single constiutent -> lambda = 0\n";
      }
      else {
        const ATOOLS::Vec4D axis(0., m_WTA ? alg->second->jetAxes(GROOM)
                                 : alg->second->jetVectors(GROOM));
      
        
        double norm = 0.;
        msg_Debugging()<<"Final States:\n";
        for(int i=0; i<ip.size(); i++) msg_Debugging()<<fl[i]<<" "
                                                      <<ip[i]<<" "
                                                      <<ip[i].PPerp()<<"\n";
        msg_Debugging()<<"... in jet:\n";
        for (const ATOOLS::Vec4D& con: constits) {
          // even with grooming , pt is defined without grooming
          const double z = con.PPerp();// / alg->second->jetVectors(0).PPerp();
          norm += z;
          const double dR = con.DR(axis) / m_R;
          
          msg_Debugging()<<con<<"\n"<<axis<<"\n";
          msg_Debugging()<<con<<" -> z = "<<con.PPerp()<< " / "<< alg->second->jetVectors(GROOM).PPerp()<<" = "<<z<<"\n";
          msg_Debugging()<<"DeltaR / R0 = "<<con.DR(axis)<<" / "<<m_R<<" = "<< dR<<" -> (DeltaR/R0)^(alpha = "<<m_alpha<<")"
                         <<" = "<<pow(dR, m_alpha)<<"\n";
          msg_Debugging()<<"z*(DeltaR/R0)^alpha = "<<z * pow(dR, m_alpha)<<"\n";
          
          lambda +=  z * pow(dR, m_alpha);
        }
        //norm = alg->second->jetVectors(0).PPerp();
        lambda /= norm;
        msg_Debugging()<<"\n";
      }
      // if(m_alpha==0.5 and GROOM==0 and lambda < 0.47 and lambda > 0.39) {
      //   if(GROOM==0) msg_Out()<<"Ungroomed\n";
      //   else msg_Out()<<"Groomed\n";
      //   msg_Out()<<"Internal: Angularity(alpha = "<<m_alpha<<") = "<<lambda<<".\n";
      //   msg_Debugging()<<"Final States:\n";
      //   for(int i=0; i<ip.size(); i++) msg_Out()<<fl[i]<<" "
      //                                           <<ip[i]<<" "
      //                                           <<ip[i].PPerp()<<"\n";
      //   msg_Out()<<lambda<<" "<<lambdaFJ<<"\n";
      //   msg_Out()<<"\n";

      // }
      return lambda;
    }           

  private:
    double m_R;
    double m_alpha;
    std::string m_pstring;
    bool m_WTA = false;
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

