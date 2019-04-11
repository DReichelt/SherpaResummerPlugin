#include "Analysis/Observable_Base.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_Limits.H"
#include "FFunction/FFunctions.H"
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
      return Obs_Params(2.0,0.0,0.0,0.0);
    }

    std::function<double(double)> FFunction(const std::vector<ATOOLS::Vec4D>& p,
                                            const std::vector<ATOOLS::Flavour>& fl) {
      if(!p_F) {
        p_F.reset(new FFUNCTION::FFunction(Name()+".dat"));
      }
      return *p_F;
    }
  
    
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

    double _Value(const std::vector<Vec4D>& ip,
                  const std::vector<Flavour>& fl,
                  const size_t &nin,
                  int njets)
    {
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
      while (n>njets) {
	if (ii!=jj) p[imap[jj]]+=p[imap[ii]];
	else THROW(fatal_error,"Invalid clustering");
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
      return dmin/Q2;
    }

    
    double _LogArg_Dynamic(double v,
                           const std::vector<Vec4T>& p,
                           const std::vector<ATOOLS::Flavour>& fl) {
      return _LogArg_Fixed(v,_Value(p,fl,2,NJETS-1));
    }

    double _ObsVal_Dynamic(double logarg,
                           const std::vector<Vec4T>& p,
                           const std::vector<ATOOLS::Flavour>& fl) {
      return _ObsVal_Fixed(logarg, _Value(p,fl,2,NJETS-1));
    }
  private:
    FFUNCTION::FFunction::Ptr p_F = nullptr;

  };// end of class YN_Durham

}// end of namespace RESUM

using namespace RESUM;

typedef YN_Durham<3> Y3_Durham;
DECLARE_GETTER(Y3_Durham,"Y3_Durham",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,Y3_Durham>::
operator()(const Parameter_Type &args) const 
{ return new Y3_Durham(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,Y3_Durham>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"Y3_Durham"; }


typedef YN_Durham<4> Y4_Durham;
DECLARE_GETTER(Y4_Durham,"Y4_Durham",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,Y4_Durham>::
operator()(const Parameter_Type &args) const 
{ return new Y4_Durham(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,Y4_Durham>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"Y4_Durham"; }


typedef YN_Durham<5> Y5_Durham;
DECLARE_GETTER(Y5_Durham,"Y5_Durham",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,Y5_Durham>::
operator()(const Parameter_Type &args) const 
{ return new Y5_Durham(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,Y5_Durham>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"Y5_Durham"; }

typedef YN_Durham<6> Y6_Durham;
DECLARE_GETTER(Y6_Durham,"Y6_Durham",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,Y6_Durham>::
operator()(const Parameter_Type &args) const 
{ return new Y6_Durham(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,Y6_Durham>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"Y6_Durham"; }
