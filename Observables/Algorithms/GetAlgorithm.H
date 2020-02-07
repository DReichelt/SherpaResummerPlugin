#include "Observables/Algorithms/Algorithm.H"
#include "Analysis/Observable_Base.H"

#include "ThrustFinder.H"

#ifdef USING_FJCONTRIB
#include "FastjetAlg.H"
#endif

namespace RESUM {

  template <typename T>
  static typename Algorithm<T>::Ptr GetAlgorithm(const std::string& name,
                                                 const std::vector<ATOOLS::Vec4<T>>& p,
                                                 const std::vector<ATOOLS::Flavour>& fl,
                                                 const size_t &nin) {
    return GetAlgorithm<T>({name},p,fl,nin);
  }
  
  template <typename T>
  static typename Algorithm<T>::Ptr GetAlgorithm(const Observable_Key& key,
                                                 const std::vector<ATOOLS::Vec4<T>>& p,
                                                 const std::vector<ATOOLS::Flavour>& fl,
                                                 const size_t &nin) {
    if(key.Name()=="ThrustFinder") {
      return std::make_shared<ThrustFinder<T>>(p,fl,nin);
    }
#ifdef USING_FJCONTRIB
    if(key.Name()=="FJmaxPTjet") {
      if(std::is_same<T,double>::value) {
        return std::make_shared<FJmaxPTjet>(p,fl,nin,key);
      }
      else {
        THROW(fatal_error, "Fastjet only available for double precision.");
        return nullptr;
      }
    }
#endif
    else {
      THROW(fatal_error,"Algorithm "+key.Name()+" not implemented!");
      return nullptr;
    }
  }
}

