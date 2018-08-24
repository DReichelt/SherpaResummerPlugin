#ifndef RESUM_THRUSTFINDER_H
#define RESUM_THRUSTFINDER_H

#include "Algorithm.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace ATOOLS;

namespace RESUM {

  template <typename T>
  inline bool mod2Cmp(const Vec3<T>& a, const Vec3<T>& b) {
    return a.Sqr() > b.Sqr();
  }

  template <typename T>
  inline T abs(const T& val) {
    return val > T(0) ? val : -val;
  }

  template <typename T>
  inline T epsilon(T ref) { return T("1e-2000"); };

  template <>
  inline double epsilon(double ref) { return 1e-12; }
  
  template <typename T>
  class ThrustFinder : public Algorithm {
  public:

    typedef std::shared_ptr<ThrustFinder> Ptr;

    ThrustFinder(const std::vector<Vec4<T>>& p,
                 const std::vector<Flavour>& fl,
                 const size_t &nin)   {
      std::vector<Vec4<T>> moms = {p.begin()+nin,p.end()};
      _calcThrust(moms,nin);
    }
  
    
    T Thrust() const { return m_t; }
    T Tau() const { return 1.-m_t; }
    const Vec3<T>& Axis() const { return m_axis; }

    
  private:

    T _calcThrust(const std::vector<Vec4<T>>& momenta, const size_t nin=2) {
      std::vector<Vec3<T>> p;
      Vec4<T> Q = {T(0),T(0),T(0),T(0)};
      for(auto&& mom: momenta) {
        p.emplace_back(mom);
        Q += mom;
      }
      if(p.size() == 2) {
        m_t = 1;
        m_axis = p[0]/p[0].Abs();
        return m_t; }
      std::vector<Vec3<T>> axes(4,{T(0),T(0),T(0)});
      std::vector<T> tval(4,T(0));
      std::sort(p.begin(), p.end(), mod2Cmp<T>);
      for(int i=0; i<4; i++) {
        int sign = i;
        for(unsigned int k=0; k<3; k++) {
          (sign % 2) == 1 ? axes.at(i) +=  p.at(k) : axes.at(i) -= p.at(k);
          sign /= 2;
        }
        axes.at(i) *= 1./axes.at(i).Abs();
        
        T diff = T(999);


        while(diff>epsilon<T>(T(0))) {
          //Todo: is this relevant for multiple precission??
          Vec3<T> foo = axes.at(i);
          axes.at(i) = {0,0,0};
          for(auto&& mom: p) {
            foo*mom > 0 ? axes.at(i) += mom : axes.at(i) -= mom;
          }
          axes.at(i) *= 1./axes.at(i).Abs();
          diff = (foo-axes.at(i)).Abs();
        }
        for(auto&& mom: p) {
          //std::cout << tval.at(i) << std::endl;
          // std::cout << std::setprecision(500) << std::endl;
          
          tval.at(i) += RESUM::abs(axes.at(i)*mom)/Q.Abs();
        }
        //std::cout << tval.at(i) << std::endl << std::endl;
      }
      m_t = 0;
      for(unsigned int i=0; i<axes.size(); i++) {
        //std::cout << std::setprecision(10) << std::endl;
        //std::cout << tval.at(i) << " " << axes.at(i) << std::endl;
        if(tval.at(i) > m_t) {
          m_t = tval.at(i);
          m_axis = axes.at(i);
        }
      }
      m_jets = std::vector<std::set<size_t>>(2);
      for(size_t i=0; i<p.size(); i++) {
        if(Vec3<T>(p[i])*m_axis < 0) m_jets.at(0).insert(i+nin);
        else m_jets.at(1).insert(i+nin);
      }
      return m_t;
    }


    
    Vec3<T> m_axis = {0,0,0};
    T m_t = 1;
  };

  
  
  
  // template <typename T>
  // class ThrustFinder : public Algorithm {
  // public:
  //   typedef std::shared_ptr<ThrustFinder> Ptr;

  //   ThrustFinder(const std::vector<Vec4<T>>& p,
  //                const std::vector<Flavour>& fl,
  //                const size_t &nin)   {
  //     size_t n = p.size();
  //     Vec3<T> lastaxis, curraxis, thrustaxis, maxaxis;
  //     T maxthrust=0., lastthrust , currthrust, thrust;
  //     std::vector<Vec3<T>> vectors(n-nin);
  //     for (size_t i(nin);i<n;++i) vectors[i-nin]=Vec3<T>(p[i][1],p[i][2],0.0);
  //     for (int pass=0;pass<2;++pass) {
  //       if (pass==1) RotateMoms(vectors,thrustaxis);
  //       sort(vectors.begin(),vectors.end(),&bigger);
  //       std::vector<Vec3<T>> initialaxes;
  //       for(unsigned int i=1;i<=(1<<(vectors.size()-1));++i) {
  //         Vec3<T> axis;
  //         for(unsigned int j=1;j<=vectors.size();++j) {
  //           int addsign=-1;
  //           if ((1<<j)*((i+(1<<(j-1))-1)/(1<<j)) >= i) addsign=1;
  //           // TODO: revise this cast
  //           axis=axis+T(addsign)*vectors[j-1];
  //         }
  //         initialaxes.push_back(axis);
  //       }
  //       sort(initialaxes.begin(),initialaxes.end(),&bigger);
  //       for(unsigned int j=0;j<initialaxes.size();j++)
  //         initialaxes[j]=initialaxes[j]/initialaxes[j].Abs();
  //       unsigned int ident=0;
  //       T sump=SumP(vectors);
  //       maxthrust=0.;
  //       for(size_t j=0;(j<initialaxes.size())&&(ident<2);j++) {
  //         curraxis=initialaxes[j];
  //         currthrust=SumNP(vectors,curraxis)/sump;
  //         lastthrust=0.0;
  //         while (currthrust>lastthrust+rpa->gen.Accu()) {
  //           lastthrust=currthrust;
  //           lastaxis=curraxis;
  //           curraxis=NewAxis(vectors,curraxis);
  //           currthrust=SumNP(vectors,curraxis)/sump;
  //         }
  //         if (lastthrust<maxthrust-rpa->gen.Accu()) break;
  //         if (lastthrust>maxthrust+rpa->gen.Accu()) {
  //           ident=0;
  //           maxthrust = lastthrust;
  //           maxaxis = lastaxis;
  //         }
  //         ident++;
  //       }
  //       if (pass==0) {
  //         m_thrust = maxthrust;
  //         m_axis = maxaxis;
  //       }
  //       m_jets = std::vector<std::set<size_t>>(2);
  //       for(size_t i=0; i<p.size(); i++) {
  //         if(Vec3<T>(p[i])*m_axis < 0) m_jets.at(0).insert(i);
  //         else m_jets.at(1).insert(i);
  //       }
  //     }
  //   }

  //   T Thrust() const { return m_thrust; }
  //   T Tau() const { return 1.-m_thrust; }
  //   const Vec3<T>& Axis() const { return m_axis; }
      
  // private:

  //   T m_thrust;
  //   Vec3<T> m_axis;

    
  //   void RotateMoms(std::vector<Vec3<T>>& ps,const Vec3<T> &ref)
  //   {
  //     for(Vec3<T>& p: ps) p -= ref*(ref*p);
  //     // for(std::vector<Vec3<T>>::iterator i(p.begin());
  //     //     i!=p.end();++i) *i=*i-ref*(ref**i);
  //   }
    
  //   Vec3<T> NewAxis(const std::vector<Vec3<T>> &p,const Vec3<T> &ref)
  //   {
  //     Vec3<T> next(0.,0.,0.);
  //     for (size_t i(0);i<p.size();++i) next+=(ref*p[i]<0.0)?-p[i]:p[i];
  //     return next/next.Abs();
  //   }
    
  //   T SumP(const std::vector<Vec3<T>> &p)
  //   { 
  //     T sum_p(0.0);
  //     for (size_t i(0);i<p.size();++i) sum_p+=p[i].Abs();
  //     return sum_p;
  //   }

  //   T SumNP(const std::vector<Vec3<T>> &p,const Vec3<T> &n)
  //   { 
  //     T sum_np(0.0);
  //     for (size_t i(0);i<p.size();++i) sum_np+=abs(p[i]*n);
  //     return sum_np;
  //   }

  //   static bool bigger(const Vec3<T> &a,const Vec3<T> &b)
  //   {
  //     return a.Sqr()>b.Sqr();
  //   }
  // };
}

#endif
