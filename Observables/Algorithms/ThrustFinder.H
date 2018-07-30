#ifndef RESUM_THRUSTFINDER_H
#define RESUM_THRUSTFINDER_H

#include "Algorithm.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace ATOOLS;

namespace RESUM {
  
  class ThrustFinder : public Algorithm {
  public:
    typedef std::shared_ptr<ThrustFinder> Ptr;
    ThrustFinder(const std::vector<Vec4D>& p,
                 const std::vector<Flavour>& fl,
                 const size_t &nin)   {
      size_t n = p.size();
      Vec3D lastaxis, curraxis, thrustaxis, maxaxis;
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
	    maxthrust = lastthrust;
            maxaxis = lastaxis;
          }
          ident++;
        }
        if (pass==0) {
          m_thrust = maxthrust;
          m_axis = maxaxis;
        }
        m_jets = std::vector<std::set<size_t>>(2);
        for(size_t i=0; i<p.size(); i++) {
          if(Vec3D(p[i])*m_axis < 0) m_jets.at(0).insert(i);
          else m_jets.at(1).insert(i);
        }
      }
    }

    double Thrust() const { return m_thrust; }
    double Tau() const { return 1.-m_thrust; }
    const Vec3D& Axis() const { return m_axis; }
      
  private:

    double m_thrust;
    Vec3D m_axis;

    
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
  };
}

#endif
