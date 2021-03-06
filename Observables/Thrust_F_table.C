#include "Analysis/Observable_Base.H"
#include "ATOOLS/Math/Poincare.H"
#include <vector>       
#include <algorithm>       


#include "FFunction/FFunctions.H"
using std::string;

using namespace ATOOLS;

namespace RESUM {

  class Thrust_F_table: public Observable_Base {
  public:

    Thrust_F_table(const Observable_Key &args): 
    Observable_Base(args) {}

    Obs_Params Parameters
      (const std::vector<ATOOLS::Vec4D>& p,
       const std::vector<ATOOLS::Flavour>& fl,
       const size_t &l) {
      return Obs_Params(1.0,1.0,0.0,0.0);
    }


    std::function<double(double,double&)> FFunction(const std::vector<Vec4D>& p,
                                            const std::vector<ATOOLS::Flavour>& fl) {
      
      if(!p_F) {
        p_F.reset(new FFUNCTION::FFunction("Additiv/Additiv_F.dat"));
      }
      return *p_F;
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

    double Value(const std::vector<Vec4D>& ip,
                 const std::vector<Flavour>& fl,
		 const size_t &nin)
    {
      size_t n = ip.size();
      // need boost to match resummation at small \tau (~1e-4)
      Vec4D sum;
      Vec4D_Vector p(&ip[0],&ip[n]);
      for (size_t i(nin);i<n;++i) sum+=p[i];
      Poincare cms(sum);
      for (size_t i(0);i<n;++i) cms.Boost(p[i]);
      // generic form as in CALT-68-836, copied from Sherpa
      Vec3D lastaxis, curraxis, thrustaxis;
      double maxthrust=0., lastthrust , currthrust, thrust;
      std::vector<Vec3D> vectors(n-nin);
      for (size_t i(nin);i<n;++i) vectors[i-nin]=Vec3D(p[i]);
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

  private:
    FFUNCTION::FFunction::Ptr p_F = nullptr;

    
  };
  

}// end of namespace RESUM

using namespace RESUM;

DECLARE_GETTER(Thrust_F_table,"Thrust_F_table",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,Thrust_F_table>::
operator()(const Parameter_Type &args) const 
{ return new Thrust_F_table(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,Thrust_F_table>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"Thrust_F_table"; }
