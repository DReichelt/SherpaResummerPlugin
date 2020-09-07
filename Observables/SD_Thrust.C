#include "Analysis/Observable_Base.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Exception.H"
#include "FFunction/FFunctions.H"
#include <algorithm>       

#ifdef USING_FJCONTRIB
#include "fastjet/contrib/SoftDrop.hh"
#endif

using namespace ATOOLS;
using namespace std;

namespace RESUM {

  class SD_Thrust: public Observable_Base {
  private:
    /// The thrust scalars.
    vector<double> m_thrusts;

    /// The thrust axes.
    vector<Vec3D> m_thrustAxes;
    
    double m_ymax;
    
  public:

    SD_Thrust(const Observable_Key &args): 
    Observable_Base(args) {
        m_zcut = to_type<double>(args.KwArg("zcut","0.0"));
        m_beta = to_type<double>(args.KwArg("beta","0"));
        m_R0 = to_type<double>(args.KwArg("R0","1"));
        m_ymax = to_type<double>(args.KwArg("ymax","-1"));
        if(m_zcut==0.) m_gmode = GROOM_MODE::NONE;
        else m_gmode = GROOM_MODE::SD;
        DEBUG_FUNC(Name()+" -> "+Tag());
        DEBUG_VAR(m_zcut);
        DEBUG_VAR(m_beta);
        DEBUG_VAR(m_R0);
        DEBUG_VAR(m_ymax);
        DEBUG_VAR(m_gmode);
    }

    Obs_Params Parameters
      (const std::vector<ATOOLS::Vec4D>& p,
       const std::vector<ATOOLS::Flavour>& fl,
       const size_t &l) {
      const double a=1.0;
      const double b=l<2?0.0:1.0;
      const double sinth=sqrt(2.0*p[2].PPerp2()/(p[0]*p[1]));
      double d=log(l<2?1.0/sinth:1.0/sqr(sinth));
      const double G_cat=0.915965594177;
      if (l<2) d+=-4.0*G_cat/M_PI-log(2.0);
      else d+=-2.0*log(2.0);
      return Obs_Params(a,b,d,0.0);
    }

    std::function<double(double, double&)> FFunction(const std::vector<ATOOLS::Vec4D>& p,
                                                     const std::vector<ATOOLS::Flavour>& fl,
                                                     const RESUM::Params& params) {
      return FFUNCTION::Additive;
    }

    double LogFac(ATOOLS::Cluster_Amplitude* ampl) {
      if(!m_unitNorm) return RESUM::Observable_Base::LogFac();
      double avg = 0;
      double n = 0;
      for(int i=(m_gmode & GROOM_MODE::SD) ? 2 : 0; i<ampl->Legs().size(); i++) {
        if(ampl->Leg(i)->Flav().StrongCharge() != 0) {
          n++;
          avg += RESUM::Observable_Base::Parameters(ampl,i).m_logdbar;
        }
      }
      avg /= n;
      return exp(avg)*RESUM::Observable_Base::LogFac();
    }

    GROOM_MODE GroomMode(double v, ATOOLS::Cluster_Amplitude* ampl,
                         const size_t &l) {
      if(l<2) {
        return m_gmode;
      }
      else {
        return GROOM_MODE::SD_COLL;
        if(v > GroomTransitionPoint(ampl, l)) return GROOM_MODE::NONE;
        else return m_gmode;        
      }
    }

//     double GroomTransitionPoint(ATOOLS::Cluster_Amplitude* ampl,
//                                 const size_t &l) {
//       
//       std::vector<Vec4T> moms(ampl->Legs().size());
//       for (size_t i=0; i<ampl->Legs().size(); ++i) {
//         moms[i]=i<ampl->NIn()?-ampl->Leg(i)->Mom():ampl->Leg(i)->Mom();
//       }
//       double sinth=sqrt(2.0*moms[2].PPerp2()/(moms[0]*moms[1]));
//       
//       Obs_Params para = RESUM::Observable_Base::Parameters(ampl,l);
//       double logfac = LogFac(ampl);
//       
//       return pow(2,m_beta)*m_zcut/pow(m_R0*sinth,m_beta)*exp(para.m_logdbar)/logfac;
//     }
    
    void RotateMoms(std::vector<Vec3D> &p,const Vec3D &ref)
    {
      for(std::vector<Vec3D>::iterator i(p.begin());
	  i!=p.end();++i) *i=*i-ref*(ref**i);
    }

    Vec3D NewAxis(const std::vector<Vec3D> &p,const Vec3D &ref)
    {
      Vec3D next(0.,0.,0.);
      for (size_t i(0);i<p.size();++i) next+=(ref*p[i]<0.0)?-p[i]:p[i];
      return next/next.Abs();
    }

    double SumP(const std::vector<Vec3D> &p)
    { 
      double sum_p(0.0);
      for (size_t i(0);i<p.size();++i) sum_p+=p[i].Abs();
      return sum_p;
    }

    double SumNP(const std::vector<Vec3D> &p,const Vec3D &n)
    { 
      double sum_np(0.0);
      for (size_t i(0);i<p.size();++i) sum_np+=dabs(p[i]*n);
      return sum_np;
    }

    static bool bigger(const Vec3D &a,const Vec3D &b)
    {
      return a.Sqr()>b.Sqr();
    }
    
    static double rapidity(const Vec4D p)
    {
      return log((p[0]+p[3])/(p[0]-p[3]))/2.;
    }

    double Value(const std::vector<Vec4D>& momenta,
                 const std::vector<Flavour>& fl,
		 const size_t &nin) {
#ifdef USING_FJCONTRIB
      size_t n = momenta.size();
      Vec3D lastaxis, curraxis, thrustaxis, maxaxis, axis;
      double maxthrust=0., lastthrust , currthrust, thrust;
      std::vector<Vec3D> vectors(n-nin);
      std::vector<Vec4D> momenta_ymax;
      momenta_ymax.push_back(momenta[0]);
      momenta_ymax.push_back(momenta[1]);
      for (size_t i(nin);i<n;++i){
          if(m_ymax<0. || abs(rapidity(momenta[i]))<m_ymax){
            vectors[i-nin]=Vec3D(momenta[i][1],momenta[i][2],0.0);
            momenta_ymax.push_back(momenta[i]);
          }
      }
      for (int pass=0;pass<2;++pass) {
	if (pass==1) RotateMoms(vectors,thrustaxis);
	sort(vectors.begin(),vectors.end(),&bigger);
	std::vector<Vec3D> initialaxes;
	for(unsigned int i=1;i<=(1<<(vectors.size()-1));++i) {
	  Vec3D axis;
	  for(unsigned int j=1;j<=vectors.size();++j) {
	    int addsign=-1;
	    if ((1<<j)*((i+(1<<(j-1))-1)/(1<<j)) >= i) addsign=1;
	    axis=axis+addsign*vectors[j-1];
	  }
	  initialaxes.push_back(axis);
	}
	sort(initialaxes.begin(),initialaxes.end(),&bigger);
	for(unsigned int j=0;j<initialaxes.size();j++)
	  initialaxes[j]=initialaxes[j]/initialaxes[j].Abs();
	unsigned int ident=0;
	double sump=SumP(vectors);
	maxthrust=0.;
	for(unsigned int j=0;(j<initialaxes.size())&&(ident<2);j++) {
	  curraxis=initialaxes[j];
	  currthrust=SumNP(vectors,curraxis)/sump;
	  lastthrust=0.0;
	  while (currthrust>lastthrust+rpa->gen.Accu()) {
	    lastthrust=currthrust;
	    lastaxis=curraxis;
	    curraxis=NewAxis(vectors,curraxis);
	    currthrust=SumNP(vectors,curraxis)/sump;
	  }
	  if (lastthrust<maxthrust-rpa->gen.Accu()) break;
	  if (lastthrust>maxthrust+rpa->gen.Accu()) {
	    ident=0;
	    maxthrust=lastthrust;
            maxaxis=lastaxis;
	  }
	  ident++;
	}
	if (pass==0) {
          thrust=maxthrust;
          axis = maxaxis;
        }
      }

      
      double ptot_perp_ungroomed = 0; //scaling factor

      
      // Partition event into hemispheres from the thrust axis
      vector<fastjet::PseudoJet> left_hemi, right_hemi;
      for(int i = nin; i < momenta_ymax.size(); i++) {
        if(!fl[i].Strong()) continue;
        Vec3D p = {momenta_ymax[i][1], momenta_ymax[i][2],0};
        ptot_perp_ungroomed += p.Abs();
        if ( axis*p >= 0) {
          right_hemi.push_back(fastjet::PseudoJet(momenta_ymax[i][1], momenta_ymax[i][2], momenta_ymax[i][3], momenta_ymax[i][0])); 
        }											
        else if ( axis*p < 0) {
          left_hemi.push_back(fastjet::PseudoJet(momenta_ymax[i][1], momenta_ymax[i][2], momenta_ymax[i][3], momenta_ymax[i][0]));
        }
        else {
          THROW(fatal_error, "Something went wrong.");
        }
      }

      if (left_hemi.size() == 0 or right_hemi.size() == 0) {
//         THROW(fatal_error,"Empty hemisphere.");
          std::cout << left_hemi.size()<< " " << right_hemi.size() << std::endl;
          std::cout << axis[0] << " " << axis[1] << std::endl;
          std::cout << momenta_ymax[0][1] << " " << momenta_ymax[0][2] << " " << momenta_ymax[0][3] << std::endl;
          std::cout << momenta_ymax[1][1] << " " << momenta_ymax[1][2] << " " << momenta_ymax[1][3] << std::endl;
          std::cout << momenta_ymax[2][1] << " " << momenta_ymax[2][2] << " " << momenta_ymax[2][3] << std::endl;
          std::cout << momenta_ymax[3][1] << " " << momenta_ymax[3][2] << " " << momenta_ymax[3][3] << std::endl;
          std::cout << momenta_ymax[4][1] << " " << momenta_ymax[4][2] << " " << momenta_ymax[4][3] << std::endl;
          return 0.;
        // empty_hemisphere += 1;
        // cout << "Empty hemisphere!! in event " << nevent << ". Veto Event and skip. This has occured: " << empty_hemisphere << " times        " << endl;
        // AddZero(ncount,0);
      }
      
      int number_of_particles = left_hemi.size() + right_hemi.size();
      
      fastjet::JetDefinition jet_def(fastjet::genkt_algorithm, 1000., 0.0, fastjet::E_scheme, fastjet::Best);
      fastjet::ClusterSequence   cs_l(left_hemi, jet_def);
      fastjet::ClusterSequence   cs_r(right_hemi, jet_def);
      vector<fastjet::PseudoJet> right_jets = cs_r.exclusive_jets(1);
      vector<fastjet::PseudoJet> left_jets = cs_l.exclusive_jets(1);
      

      fastjet::contrib::SoftDrop softdrop(m_beta, m_zcut, m_R0); // Inputs are beta, zcut, R0
      vector<fastjet::PseudoJet> sd_right_hemi = softdrop(right_jets);
      vector<fastjet::PseudoJet> sd_left_hemi = softdrop(left_jets);
      vector<fastjet::PseudoJet>  right = sd_right_hemi[0].constituents();
      vector<fastjet::PseudoJet>  left = sd_left_hemi[0].constituents();
      vector<fastjet::PseudoJet>  left_ungroom = left_jets[0].constituents();
      vector<fastjet::PseudoJet>  right_ungroom = right_jets[0].constituents();

      // Calculate ptot
      double ptot = 0;
      
      for (unsigned i = 0; i < right.size(); i++){
        ptot += right[i].pt();
      }

      for (unsigned i = 0; i < left.size(); i++){
        ptot += left[i].pt();
      } 
      
      // Calculate SD thrust axes
      Vec3D p_r = {sd_right_hemi[0].px(), sd_right_hemi[0].py(), 0};
      Vec3D p_l = {sd_left_hemi[0].px(), sd_left_hemi[0].py(), 0};
      
      // thrust
      double T = ( p_l.Abs() + p_r.Abs() );
      double tau = (ptot - T) / ptot_perp_ungroomed;
      return tau;  
#else
      THROW(fatal_error, "No grooming without fastjet.");
      return 0;
#endif
    }
    
  };// end of class Thrust

}// end of namespace RESUM

using namespace RESUM;

DECLARE_GETTER(SD_Thrust,"SD_Thrust",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,SD_Thrust>::
operator()(const Parameter_Type &args) const 
{ return new SD_Thrust(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,SD_Thrust>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"SD_Thrust"; }
