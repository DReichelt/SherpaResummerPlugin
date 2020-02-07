#ifndef RESUM_ALGORITHM_H
#define RESUM_ALGORITHM_H
#include <memory>
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Vector.H"

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
    inline const std::vector<ATOOLS::Vec4<T>>&  jetVectors() const {return m_jetVectors;}
    inline ATOOLS::Vec4<T> jetVectors(size_t i) const {return m_jetVectors.at(i);}
    inline const std::vector<T> jetMasses() const {return m_jetMasses;}
    inline T jetMasses(size_t i) const {return m_jetMasses.at(i);}
    inline const std::vector<T> jetBroads() const {return m_jetBroads;}
    inline T jetBroads(size_t i) const {return m_jetBroads.at(i);}
    inline const std::vector<T> jetScales() const {return m_jetScales;}
    inline T jetScales(size_t i) const {return m_jetScales.at(i);}
    inline T Q2() const {return m_Q2;}
    
  protected:
    std::vector<std::set<size_t>> m_jets;
    std::vector<ATOOLS::Vec3<T>> m_jetAxes;
    std::vector<ATOOLS::Vec4<T>> m_jetVectors;
    std::vector<T> m_jetMasses;
    std::vector<T> m_jetBroads;
    std::vector<T> m_jetScales;
    T m_Q2;
  };
}


#endif
