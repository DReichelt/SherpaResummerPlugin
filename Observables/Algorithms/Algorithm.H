#ifndef RESUM_ALGORITHM_H
#define RESUM_ALGORITHM_H
#include <memory>
#include "ATOOLS/Org/Exception.H"

namespace RESUM {
    class Algorithm {
    public:
      typedef std::shared_ptr<Algorithm> Ptr;
      Algorithm() {}

      const std::vector<std::set<size_t>>& Jets() const { return m_jets; }
      
    protected:
      std::vector<std::set<size_t>> m_jets;
    };
}


#include "ThrustFinder.H"


namespace RESUM {
  static Algorithm::Ptr GetAlgorithm(const std::string& name,
                                     const std::vector<ATOOLS::Vec4D>& p,
                                     const std::vector<ATOOLS::Flavour>& fl,
                                     const size_t &nin) {
    if(name=="ThrustFinder") return std::make_shared<ThrustFinder>(p,fl,nin);
    else {
      THROW(fatal_error,"Algorithm "+name+" not implemented!");
      return nullptr;
    }
  }
}


#endif
