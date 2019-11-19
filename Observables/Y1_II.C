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

  class Y1_II: public Observable_Base {
  public:

    Y1_II(const Observable_Key &args): 
    Observable_Base(args) {}

    Obs_Params Parameters
      (const std::vector<ATOOLS::Vec4D>& p,
       const std::vector<ATOOLS::Flavour>& fl,
       const size_t& l) {
      return Obs_Params(2.0,0.0,0.0,0.0);
    }

    std::function<double(double,double&)> FFunction(const std::vector<ATOOLS::Vec4D>& p,
                                                    const std::vector<ATOOLS::Flavour>& fl, 
                                                    const RESUM::Params& params) {
      return FFUNCTION::Additive;
    }


    
    double KT2(const Vec4D &p1, const Vec4D &p2=Vec4D()) const
    {
      if (p2==Vec4D()) return p1.PPerp2();
      return 2.0*sqr(Min(p1.PPerp2(),p2.PPerp2()))*
	(cosh(p1.Y()-p2.Y())-cos(p1.DPhi(p2)));
    }

    double Value(const std::vector<Vec4D>& ip,
                 const std::vector<Flavour>& fl,
		 const size_t &nin)
    {
      size_t nn = ip.size();
      Vec4D Q;
      Vec4D_Vector p;
      for (size_t i(nin);i<nn;++i) {
	if (fl[i].Strong()) p.push_back(ip[i]);
	else Q+=ip[i];
      }
      std::vector<int> imap(p.size());
      for (int i=0;i<imap.size();++i) imap[i]=i;
      std::vector<std::vector<double> > kt2ij
	(imap.size(),std::vector<double>(imap.size()));
      // calc first matrix
      int ii=0, jj=0, n=p.size();
      double dmin=KT2(p[0]);
      for (int i=0;i<n;++i) {
	double di=kt2ij[i][i]=KT2(p[i]);
	if (di<dmin) { dmin=di; ii=i; jj=i; }
	for (int j=0;j<i;++j) {
	  double dij=kt2ij[i][j]=KT2(p[i],p[j]);
	  if (dij<dmin) { dmin=dij; ii=i; jj=j; }
	}
      }
      // recalc matrix
      while (n>1) {
	int jjx=imap[jj];
	if (ii!=jj) { p[jjx]+=p[imap[ii]]; p[jjx][0]=p[jjx].PSpat(); }
	--n;
	for (int i=ii;i<n;++i) imap[i]=imap[i+1];
	kt2ij[jjx][jjx]=KT2(p[jjx]);
	for (int j=0;j<jj;++j) kt2ij[jjx][imap[j]]=KT2(p[jjx],p[imap[j]]);
	for (int i=jj+1;i<n;++i) kt2ij[imap[i]][jjx]=KT2(p[imap[i]],p[jjx]);
	ii=jj=0;
	dmin=kt2ij[imap[0]][imap[0]];
	for (int i=0;i<n;++i) {
	  int ix=imap[i];
	  double di=kt2ij[ix][ix];
	  if (di<dmin) { dmin=di; ii=jj=i; }
	  for (int j=0;j<i;++j) {
	    int jx=imap[j];
	    double dij=kt2ij[ix][jx];
	    if (dij<dmin) { dmin=dij; ii=i; jj=j; }
	  }
	}
      }
      return dmin/Q.Abs2();
    }

  };// end of class Y1_II

}// end of namespace RESUM

using namespace RESUM;

DECLARE_GETTER(Y1_II,"Y1_II",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,Y1_II>::
operator()(const Parameter_Type &args) const 
{ return new Y1_II(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,Y1_II>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"Y1_II"; }
