#include "Analysis/Observable_Base.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Poincare.H"
#include "FFunction/FFunctions.H"
#include <algorithm>       

using namespace ATOOLS;
using namespace std;

namespace RESUM {

  class Thrust: public Observable_Base {
  private:
    /// The thrust scalars.
    vector<double> m_thrusts;

    /// The thrust axes.
    vector<Vec3D> m_thrustAxes;
    
  public:

    Thrust(const Observable_Key &args): 
    Observable_Base(args) {}

    Obs_Params Parameters
      (const std::vector<ATOOLS::Vec4D>& p,
       const std::vector<ATOOLS::Flavour>& fl,
       const size_t &l) {
      double a=1.0;
      double b=l<2?0.0:1.0;
      double sinth=sqrt(2.0*p[2].PPerp2()/(p[0]*p[1]));
      double d=log(l<2?1.0/sinth:1.0/sqr(sinth));
      const double G_cat=0.915965594177;
      if (l<2) d+=-4.0*G_cat/M_PI-log(2.0);
      else d+=-2.0*log(2.0);
      return Obs_Params(a,b,d,0.0);
    }

    std::function<double(double,double&)> FFunction(const std::vector<ATOOLS::Vec4D>& p,
                                                    const std::vector<ATOOLS::Flavour>& fl, 
                                                    const RESUM::Params& params) {
      return FFUNCTION::Additive;
    }

    
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

    double Value(const std::vector<Vec4D>& p,
                 const std::vector<Flavour>& fl,
		 const size_t &nin)
    {
      size_t n = p.size();
      Vec3D lastaxis, curraxis, thrustaxis;
      double maxthrust=0., lastthrust , currthrust, thrust;
      std::vector<Vec3D> vectors(n-nin);
      for (size_t i(nin);i<n;++i) vectors[i-nin]=Vec3D(p[i][1],p[i][2],0.0);
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
	  }
	  ident++;
	}
	if (pass==0) thrust=maxthrust; 
      }
      return 1.0-thrust;
    }
    
  };// end of class Thrust

}// end of namespace RESUM

using namespace RESUM;

DECLARE_GETTER(Thrust,"Thrust",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,Thrust>::
operator()(const Parameter_Type &args) const 
{ return new Thrust(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,Thrust>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"Thrust"; }
