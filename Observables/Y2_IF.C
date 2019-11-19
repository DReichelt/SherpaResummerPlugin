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

  class Y2_IF: public Observable_Base {
  public:

    Y2_IF(const Observable_Key &args): 
    Observable_Base(args) {}

    Obs_Params Parameters
      (const std::vector<Vec4D>& p,
       const std::vector<Flavour>& fl,
       const size_t& l) {
      return Obs_Params(2.0,0.0,0.0,0.0);
    }

    std::function<double(double,double&)> FFunction(const std::vector<ATOOLS::Vec4D>& p,
                                            const std::vector<ATOOLS::Flavour>& fl) {
      // @TODO is this correct??
      return FFUNCTION::Additive;
    }

    
    double KT2(const Vec4D &p1, const Vec4D &p2) const
    {
      if (p2.PPerp()==0.0) return 2.0*sqr(p1[0])*(1.0-p1.CosTheta(p2));
      return 2.0*sqr(Min(p1[0],p2[0]))*(1.0-p1.CosTheta(p2));
    }

    double Value(const std::vector<Vec4D>& ip,
                 const std::vector<Flavour>& fl,
		 const size_t &nin)
    {
      size_t nn = ip.size();
      // identify particles, reconstruct beam
      Vec4D asum;
      Vec4D_Vector moms(1);
      for (size_t i(nin);i<nn;++i) {
	if (fl[i].Strong()) moms.push_back(ip[i]);
	asum+=ip[i];
      }
      int beam(fl[0].IsLepton());
      double pm(asum.Mass()*exp(beam?-asum.Y():asum.Y()));
      moms[0]=-Vec4D(pm,0.0,0.0,beam?-pm:pm)/2.0;
      // boost into Breit frame
      Vec4D qq, pp(rpa->gen.PBeam(beam));
      pp[0]=pp.PSpat();
      for (size_t i(0);i<moms.size();++i) qq+=moms[i];
      Poincare cms(pp+qq);
      double Q2(-qq.Abs2()), x(Min(Q2/(2.0*pp*qq),1.0));
      msg_Debugging()<<"Q^2 = "<<Q2<<"\n";
      double E(sqrt(Q2)/((beam?2.0:-2.0)*x));
      Vec4D P(sqrt(E*E+pp.Abs2()),0.0,0.0,-E);
      Vec4D q(0.0,0.0,0.0,2.0*x*E);
      cms.Boost(pp);
      cms.Boost(qq);
      Poincare zrot(pp,beam?-Vec4D::ZVEC:Vec4D::ZVEC);
      zrot.Rotate(pp);
      zrot.Rotate(qq);
      Poincare breit(P+q);
      breit.BoostBack(pp);
      breit.BoostBack(qq);
      double Q(sqrt(Q2));
      moms[0]=-moms[0];
      for (int i(0);i<moms.size();++i) {
	msg_Debugging()<<"p["<<i<<"] = "<<moms[i];
	cms.Boost(moms[i]);
	zrot.Rotate(moms[i]);
	breit.BoostBack(moms[i]);
	msg_Debugging()<<" -> "<<moms[i]<<"\n";
      }
      Vec4D_Vector p(&moms[1],&moms[moms.size()]);
      std::vector<int> imap(p.size());
      for (int i=0;i<imap.size();++i) imap[i]=i;
      std::vector<std::vector<double> > kt2ij
	(imap.size(),std::vector<double>(imap.size()));
      // calc first matrix
      int ii=0, jj=0, n=p.size();
      double dmin=KT2(p[0],P);
      for (int i=0;i<n;++i) {
	double di=kt2ij[i][i]=KT2(p[i],P);
	if (di<dmin) { dmin=di; ii=i; jj=i; }
	for (int j=0;j<i;++j) {
	  double dij=kt2ij[i][j]=KT2(p[i],p[j]);
	  if (dij<dmin) { dmin=dij; ii=i; jj=j; }
	}
      }
      // recalc matrix
      while (n>2) {
	int jjx=imap[jj];
	if (ii!=jj) { p[jjx]+=p[imap[ii]]; p[jjx][0]=p[jjx].PSpat(); }
	--n;
	for (int i=ii;i<n;++i) imap[i]=imap[i+1];
	kt2ij[jjx][jjx]=KT2(p[jjx],P);
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
      return dmin/Q2;
    }

  };// end of class Y2_IF

}// end of namespace RESUM

using namespace RESUM;

DECLARE_GETTER(Y2_IF,"Y2_IF",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,Y2_IF>::
operator()(const Parameter_Type &args) const 
{ return new Y2_IF(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,Y2_IF>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"Y2_IF"; }
