#include "Analysis/Observable_Base.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Org/Message.H"
#include "FFunction/FFunctions.H"
#include "Main/Params.H"
#include <vector>       
#include <algorithm>       

#define s_ymax std::numeric_limits<double>::max()

using namespace ATOOLS;

namespace RESUM {

  template <int NJETS>
  class YN_Durham: public Observable_Base {
  public:

    YN_Durham(const Observable_Key &args): 
    Observable_Base(args) {}


    Obs_Params Parameters
      (const std::vector<Vec4D>& p,
       const std::vector<Flavour>& fl,
       const size_t &l) {
      return {2.,0.,0.,1.};
    }


    std::function<double(double,double&)> FFunction(const std::vector<ATOOLS::Vec4D>& p,
                                                    const std::vector<ATOOLS::Flavour>& fl,
                                                    const RESUM::Params& params) {
      
      if(NJETS > 6) {
        return [](double Rp,double& exp) {exp=0.;return 1.;};
      }
      size_t N_gluon = 0;

      double NC = params.NC();
      double CF = params.CF();
      double CA = params.CA();

      double num = 0;
      double den = 0;
      for (size_t i=2; i<fl.size(); i++){
	if (fl[i].IsGluon()){
	  N_gluon += 1;
          num += CA*CA;
          den += CA;
	}
        else if(fl[i].IsQuark()) {
          num += CF*CF;
          den += CF;
        }
      }
      double F2 = -pow(M_PI,2)/8. * 0.5*num/den/den;
      if(!p_F1 and (N_gluon == 0 or N_gluon == 1)) {
        std::string name = Name() + "_" + std::to_string(N_gluon) + "g";
        if(params.largeNC()) name += "_largeNc";
        name += ".dat";
        p_F1.reset(new FFUNCTION::FFunction(name, F2));
      }
            
      if(!p_F2 and (N_gluon == 2 or N_gluon == 3)) {
        std::string name = Name() + "_" + std::to_string(N_gluon) + "g";
        if(params.largeNC()) name += "_largeNc";
        name += ".dat";
        p_F2.reset(new FFUNCTION::FFunction(name, F2));
      }
      

      if (N_gluon < 2 && p_F1 != nullptr) return *p_F1;
      if (N_gluon < 4 && p_F2 != nullptr) return *p_F2;
      THROW(fatal_error,"No F function for events with " + std::to_string(N_gluon) + " gluons.");
      return 0;
    }

    virtual std::set<size_t> ResumMult() {return {NJETS-1};}
    virtual size_t ResumQCDorderBorn() {return NJETS-3;};
    virtual size_t ResumQCDorderLO() {return NJETS-2;}
    virtual size_t ResumQCDorderNLO() {return NJETS-1;}
  
    
    double KT2(const Vec4D &p1, const Vec4D &p2) const
    {
      return 2.0*sqr(Min(p1[0],p2[0]))*(1.0-p1.CosTheta(p2));
    }

    double Value(const std::vector<Vec4D>& ip,
                 const std::vector<Flavour>& fl,
		 const size_t &nin)
    {
      return _Value(ip,fl,nin,NJETS);
    }


    int _rejected = 0;
    
    double _Value(const std::vector<Vec4D>& ip,
                  const std::vector<Flavour>& fl,
                  const size_t &nin,
                  int njets, int nmax=9999)
    {
      for(const Vec4D& p: ip) {
        if(p.Nan()) {
          _rejected++;
          msg_Error()<<METHOD<<": Nan in momenta, rejecting (returning 0)\n";
          for(const Vec4D& p: ip) msg_Error()<<p<<"\n";
          msg_Error()<<"Rejected "<<_rejected<<" events so far.\n";
          return 0;
        }
      }
      Vec4D sum;
      size_t nn = ip.size();
      Vec4D_Vector p(&ip[nin],&ip[nn]);
      if(p.size() < njets) return 0;
      for (size_t i(0);i<p.size();++i) sum+=p[i];
      Poincare cms(sum);
      for (size_t i(0);i<p.size();++i) cms.Boost(p[i]);
      double Q2(sum.Abs2());
      std::vector<int> imap(p.size());
      for (int i=0;i<imap.size();++i) imap[i]=i;
      std::vector<std::vector<double> > kt2ij
	(imap.size(),std::vector<double>(imap.size()));
      int ii=0, jj=0, n=p.size();
      double dmin=Q2;
      for (int i=0;i<n;++i)
	for (int j=0;j<i;++j) {
	  double dij=kt2ij[i][j]=KT2(p[i],p[j]);
	  if (dij<dmin) { dmin=dij; ii=i; jj=j; }
	}
      m_yij.clear();
      while (n>njets) {
	if (ii!=jj) p[imap[jj]]+=p[imap[ii]];
	else {
          msg_Error()<<"Something went wrong clustering the following: ";
          for(const Vec4D& p: ip) msg_Error()<<p<<"\n";
          THROW(fatal_error,"Invalid clustering");
        }
        if(n <= nmax) m_yij.push_back(dmin/Q2);
	--n;
	for (int i=ii;i<n;++i) imap[i]=imap[i+1];
	int jjx=imap[jj];
	for (int j=0;j<jj;++j) kt2ij[jjx][imap[j]]=KT2(p[jjx],p[imap[j]]);
	for (int i=jj+1;i<n;++i) kt2ij[imap[i]][jjx]=KT2(p[imap[i]],p[jjx]);
	ii=jj=0; dmin=Q2;
	for (int i=0;i<n;++i)
	  for (int j=0;j<i;++j) {
	    double dij=kt2ij[imap[i]][imap[j]];
	    if (dij<dmin) { dmin=dij; ii=i; jj=j; }
	  }
      }
      if(n <= nmax) m_yij.push_back(dmin/Q2);
      return dmin/Q2;
    }

    const std::vector<double>&  AllValues(const std::vector<Vec4D>& ip,
                                      const std::vector<Flavour>& fl,
                                          const size_t &nin, int nmax) {
      _Value(ip,fl,nin,NJETS,nmax);
      return m_yij;
    }
    
    double _Endpoint_Dynamic(ATOOLS::Cluster_Amplitude* ampl) {
      Vec4D_Vector moms(ampl->Legs().size());
      Flavour_Vector flavs(ampl->Legs().size());
      for (size_t i(0);i<ampl->Legs().size();++i) {
        moms[i]=i<ampl->NIn()?-ampl->Leg(i)->Mom():ampl->Leg(i)->Mom();
        flavs[i]=i<ampl->NIn()?ampl->Leg(i)->Flav().Bar():ampl->Leg(i)->Flav();
      }
      double ret = _Value(moms,flavs,ampl->NIn(),NJETS-1);
      // if (ret > 2/pow(NJETS,2)) msg_Out()<<ret<<" "<<2/pow(NJETS,2)<<"\n";
      return ret < 2/pow(NJETS,2) ? ret :  2/pow(NJETS,2);     
    }


  private:
    FFUNCTION::FFunction::Ptr p_F1 = nullptr;
    FFUNCTION::FFunction::Ptr p_F2 = nullptr;


    std::vector<double> m_yij;
  };// end of class YN_Durham

}// end of namespace RESUM
