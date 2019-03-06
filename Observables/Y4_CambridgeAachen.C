#include "Analysis/Observable_Base.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_Limits.H"
#include <vector>       
#include <algorithm>       

#define s_ymax std::numeric_limits<double>::max()

using namespace ATOOLS;

namespace RESUM {

  class Y4_CambridgeAachen: public Observable_Base {
  public:

    Y4_CambridgeAachen(const Observable_Key &args): 
    Observable_Base(args) {}

    Obs_Params Parameters
      (const std::vector<Vec4D>& p,
       const std::vector<Flavour>& fl,
       const size_t &l) {
      return Obs_Params(2.0,0.0,0.0,0.0);
    }

    double CalcF(const double Rp) {
      return 1.;
    }

    double KT2(const Vec4D &p1, const Vec4D &p2) const
    {
      return 2.0*sqr(Min(p1[0],p2[0]))*(1.0-p1.CosTheta(p2));
    }

    double y12(const Vec4D &p1, const Vec4D &p2) const {
      return (1.0-p1.CosTheta(p2));
    }
    
    double Value(const std::vector<Vec4D>& ip,
                 const std::vector<Flavour>& fl,
		 const size_t &nin)
    {
      Vec4D sum;
      size_t nn = ip.size();
      Vec4D_Vector p(&ip[nin],&ip[nn]);
      for (size_t i(0);i<p.size();++i) sum+=p[i];
      Poincare cms(sum);
      for (size_t i(0);i<p.size();++i) cms.Boost(p[i]);
      double Q2(sum.Abs2());
      std::vector<int> imap(p.size());
      for (int i=0;i<imap.size();++i) imap[i]=i;
      std::vector<std::vector<double> > yij
	(imap.size(),std::vector<double>(imap.size()));
      int ii=0, jj=0, n=p.size();
      double dmin=Q2;
      for (int i=0;i<n;++i)
	for (int j=0;j<i;++j) {
	  double dij=yij[i][j]=y12(p[i],p[j]);
	  if (dij<dmin) { dmin=dij; ii=i; jj=j; }
	}
      double kt2_4 = -1;
      while (n>4) {
	if (ii!=jj) {
          p[imap[jj]]+=p[imap[ii]];
        }
	else THROW(fatal_error,"Invalid clustering");
	--n;
	for (int i=ii;i<n;++i) imap[i]=imap[i+1];
	int jjx=imap[jj];
	for (int j=0;j<jj;++j) yij[jjx][imap[j]]=y12(p[jjx],p[imap[j]]);
	for (int i=jj+1;i<n;++i) yij[imap[i]][jjx]=y12(p[imap[i]],p[jjx]);
	ii=jj=0; dmin=Q2;
	for (int i=0;i<n;++i)
	  for (int j=0;j<i;++j) {
	    double dij=yij[imap[i]][imap[j]];
	    if (dij<dmin) {
              dmin=dij;
              ii=i;
              jj=j;
            }
	  }
      }
      kt2_4 = KT2(p[imap[jj]],p[imap[ii]]);
      return kt2_4/Q2;
    }

  };// end of class Y4_CambridgeAachen

}// end of namespace RESUM

using namespace RESUM;

DECLARE_GETTER(Y4_CambridgeAachen,"Y4_CambridgeAachen",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,Y4_CambridgeAachen>::
operator()(const Parameter_Type &args) const 
{ return new Y4_CambridgeAachen(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,Y4_CambridgeAachen>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"Y4_CambridgeAachen"; }
