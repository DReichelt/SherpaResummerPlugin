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

  class Thrust_IF: public Observable_Base {
  public:

    Thrust_IF(const Observable_Key &args): 
    Observable_Base(args) {}

    Obs_Params Parameters
    (const ATOOLS::Vec4D *p,const ATOOLS::Flavour *fl,
     const size_t &n,const size_t &l) {
      return Obs_Params(1.0,1.0,0.0,0.0);
    }

    double Value(const Vec4D *ip,const Flavour *fl,
		 const size_t &nn,const size_t &nin)
    {
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
      Vec4D p(sqrt(E*E+pp.Abs2()),0.0,0.0,-E);
      Vec4D q(0.0,0.0,0.0,2.0*x*E);
      cms.Boost(pp);
      cms.Boost(qq);
      Poincare zrot(pp,beam?-Vec4D::ZVEC:Vec4D::ZVEC);
      zrot.Rotate(pp);
      zrot.Rotate(qq);
      Poincare breit(p+q);
      breit.BoostBack(pp);
      breit.BoostBack(qq);
      Vec4D n(pp/pp.PSpat()), sum;
      Vec4D nb(Vec4D(n[0],-n[1],-n[2],-n[3]));
      double tau(0.0), Q(sqrt(Q2));
      msg_Debugging()<<"n  = "<<n<<", nb = "<<nb<<"\n";
      moms[0]=-moms[0];
      // hep-ph/9912488, Eq.(2.6)
      for (int i(0);i<moms.size();++i) {
	msg_Debugging()<<"p["<<i<<"] = "<<moms[i];
	cms.Boost(moms[i]);
	zrot.Rotate(moms[i]);
	breit.BoostBack(moms[i]);
	msg_Debugging()<<" -> "<<moms[i]<<"\n";
	sum+=i==0?-moms[i]:moms[i];
	double tb(n*moms[i]), tj(nb*moms[i]);
	if (i>0 && tb>tj) tau+=dabs(moms[i][3]);
      }
      msg_Debugging()<<"mom sum = "<<sum<<"\n";
      tau=1.0-2.0/Q*tau;
      msg_Debugging()<<"\\tau = "<<tau<<"\n";
      return tau;
    }

  };// end of class Thrust_IF

}// end of namespace RESUM

using namespace RESUM;

DECLARE_GETTER(Thrust_IF,"Thrust_IF",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,Thrust_IF>::
operator()(const Parameter_Type &args) const 
{ return new Thrust_IF(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,Thrust_IF>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"Thrust_IF"; }
