#ifndef RESUM_ALGORITHM_H
#define RESUM_ALGORITHM_H
#include <memory>
#include "ATOOLS/Org/Exception.H"

namespace RESUM {
    class Algorithm {
    public:
      typedef std::shared_ptr<Algorithm> Ptr;
      Algorithm() {}

      virtual ~Algorithm() = default;
      
      const std::vector<std::set<size_t>>& Jets() const { return m_jets; }
      template <typename T>
      std::vector<T> apply(std::vector<T> all, size_t jetIdx) {
        std::vector<T> ret(m_jets.at(jetIdx).size());
        auto idit = m_jets[jetIdx].begin();
        for(auto retit = ret.begin(); retit != ret.end(); retit++, idit++) {
          *retit = all.at(*idit);
        }
        return ret;
      }
      
    protected:
      std::vector<std::set<size_t>> m_jets;
    };
}


#include "ThrustFinder.H"


namespace RESUM {
  template <typename T>
  static Algorithm::Ptr GetAlgorithm(const std::string& name,
                                     const std::vector<ATOOLS::Vec4<T>>& p,
                                     const std::vector<ATOOLS::Flavour>& fl,
                                     const size_t &nin) {
    if(name=="ThrustFinder") return std::make_shared<ThrustFinder<T>>(p,fl,nin);
    else {
      THROW(fatal_error,"Algorithm "+name+" not implemented!");
      return nullptr;
    }
  }
}


#endif
