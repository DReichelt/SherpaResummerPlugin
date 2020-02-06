#ifndef RESUM_ALGORITHM_H
#define RESUM_ALGORITHM_H
#include <memory>
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/MathTools.H"


namespace RESUM {

  template <typename T>
  class Algorithm {
  public:
    typedef std::shared_ptr<Algorithm> Ptr;
    Algorithm() {}
    
    virtual ~Algorithm() = default;
    
    const std::vector<std::set<size_t>>& Jets() const { return m_jets; }
    
    template <typename X>
    std::vector<X> apply(std::vector<X> all, size_t jetIdx) {
      std::vector<X> ret(m_jets.at(jetIdx).size());
      auto idit = m_jets[jetIdx].begin();
      for(auto retit = ret.begin(); retit != ret.end(); retit++, idit++) {
        *retit = all.at(*idit);
      }
      return ret;
    }
    
    inline size_t numJets() const {return m_jets.size();}
    inline const std::vector<ATOOLS::Vec3<T>>&  jetAxes() const {return m_jetAxes;}
    inline ATOOLS::Vec3<T> jetAxes(size_t i) const {return m_jetAxes.at(i);}
    inline const std::vector<T> jetMasses() const {return m_jetMasses;}
    inline T jetMasses(size_t i) const {return m_jetMasses.at(i);}
    inline const std::vector<T> jetBroads() const {return m_jetBroads;}
    inline T jetBroads(size_t i) const {return m_jetBroads.at(i);}
    inline T Q2() const {return m_Q2;}
    
  protected:
    std::vector<std::set<size_t>> m_jets;
    std::vector<ATOOLS::Vec3<T>> m_jetAxes;
    std::vector<T> m_jetMasses;
    std::vector<T> m_jetBroads;
    T m_Q2;
  };
}


#include "ThrustFinder.H"
#include "FastjetAlg.H"

namespace RESUM {
  template <typename T>
  static typename Algorithm<T>::Ptr GetAlgorithm(const std::string& name,
                                                 const std::vector<ATOOLS::Vec4<T>>& p,
                                                 const std::vector<ATOOLS::Flavour>& fl,
                                                 const size_t &nin) {
    if(name=="ThrustFinder") return std::make_shared<ThrustFinder<T>>(p,fl,nin);
    // if(name=="FJmaxPTjet") {
    //   if(std::is_same<T,double>::value) {
    //     return std::make_shared<FJmaxPTjet>(p,fl,nin);
    //   }
    //   else {
    //     THROW(fatal_error, "Fastjet only available for double precision.");
    //     return nullptr;
    //   }
    // }
    else {
      THROW(fatal_error,"Algorithm "+name+" not implemented!");
      return nullptr;
    }
  }
}


#endif
