#ifndef RESUM_EPARAMETER_H
#define RESUM_EPARAMETER_H

#include "Analysis/Observable_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "Observables/Algorithms/GetAlgorithm.H"
#include "FFunction/FFunctions.H"
#include "Math/Matrix.H"

namespace RESUM {
  template <typename T>
  class EParameter_Template: public Observable_Base_Template<T> {
  private:
    using Vec4T = ATOOLS::Vec4<T>;
    using Vec3T = ATOOLS::Vec3<T>;
    
  public:
    EParameter_Template(const Observable_Key& args) : Observable_Base_Template<T>(args) {}

    Obs_Params Parameters(const std::vector<Vec4T>& moms,
                          const std::vector<ATOOLS::Flavour>& flavs,
                          const size_t& l=0) {
      return {1,0,log(2.),1};
    }

    T Value(const std::vector<Vec4T>& moms,
                 const std::vector<ATOOLS::Flavour>& flavs,
                 const size_t& nin) {
      std::map<std::string, typename Algorithm<T>::Ptr> dummy;
      return Value(moms, flavs, dummy, nin);
    }

    std::function<double(double,double&)> FFunction(const std::vector<ATOOLS::Vec4D>& p,
                                                    const std::vector<ATOOLS::Flavour>& fl, 
                                                    const RESUM::Params& params) {
      /// @TODO insert correct function
      return FFUNCTION::Additive;
    }
       
    T Value(const std::vector<Vec4T>& moms,
            const std::vector<ATOOLS::Flavour>& flavs,
            std::map<std::string, typename Algorithm<T>::Ptr>& algorithms,
            const size_t& nin) {
      if(moms.size()-nin == 2) return T(0);
      const std::string name = "ThrustFinder";
      auto alg = algorithms.find(name);
      if(alg==algorithms.end()) {
        alg = algorithms.insert({name,
                                 GetAlgorithm<T>(name, moms, flavs, nin)}).first;
        msg_Debugging()<<"Found jets.\n";
      }
      else msg_Debugging()<<"Reusing jets found earlier.\n";
      Vec3T axis = std::dynamic_pointer_cast<ThrustFinder<T>>(alg->second)->Axis();
      Vec3T n1 = {T(0),T(0),T(0)};
      std::vector<Vec3T> projMoms(moms.size());
      for(size_t i=0; i<projMoms.size()-nin; i++) {
        const Vec3T mom = Vec3T(moms.at(i+nin));
        projMoms[i] = mom-(mom*axis) * axis;
        if(projMoms[i].Abs() > n1.Abs()) n1 = projMoms[i];
      }
      n1 *= T(1)/n1.Abs();
      Vec3T n2 = cross(n1,axis);
      n2 *= T(1)/n2.Abs();
      RESUM::Matrix<T> L(2,2,0);
      for(const Vec3T& m: projMoms) {
        if(!ATOOLS::IsZero(m.Abs())) {
          L(0,0) += (m*n1) * (n1*m) / m.Abs();
          L(0,1) += (m*n1) * (n2*m) / m.Abs();
          L(1,0) += (m*n2) * (n1*m) / m.Abs();
          L(1,1) += (m*n2) * (n2*m) / m.Abs();
        }
      }
      L /= sqrt(alg->second->Q2());
      const T l1 = (L(0,0)+L(1,1) + sqrt(pow(L(0,0)+L(1,1),2)-T(4)*(L(0,0)*L(1,1)-L(0,1)*L(1,0))))/T(2);
      const T l2 = (L(0,0)+L(1,1) - sqrt(pow(L(0,0)+L(1,1),2)-T(4)*(L(0,0)*L(1,1)-L(0,1)*L(1,0))))/T(2);
      return l1 > l2 ? l1 : l2;
    }

  };

}


#endif
