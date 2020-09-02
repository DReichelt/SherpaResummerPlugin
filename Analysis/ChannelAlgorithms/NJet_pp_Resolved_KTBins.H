#ifndef RESUM_NJET_PP_RES_H
#define RESUM_NJET_PP_RES_H
#include "Analysis/ChannelAlgorithms/ChannelAlgorithm_Base.H"
#include <string>
#include <vector>
#include "ATOOLS/Math/Vector.H"


namespace RESUM {
  class NJet_pp_Resolved_KTBins : public ChannelAlgorithm_Base {
  public:
    NJet_pp_Resolved_KTBins(const ChAlg_Key& parameters);
    std::string Channel(const std::vector<ATOOLS::Vec4D>& ip,
                        const std::vector<ATOOLS::Flavour>& fl,
                        const size_t &nin,
                        std::vector<ATOOLS::Vec4D>* pout) override;

    std::vector<std::string> ChannelNames() override {return m_channelNames;}

    enum MODE {
               ALL = 0,
               BLAND = 1   << 0,
               ZPROD = 1   << 1,
               BLAND_Z = BLAND|ZPROD
    };
    
  private:
        
    MODE m_mode = MODE::ALL;
    std::vector<std::string> m_channelNames;
    std::vector<double> m_binEdges;
    std::vector<std::string> m_edges;
    std::vector<ChannelAlgorithm_Base::Ptr> m_resolvers;

  }; // end of class NJet_pp_Resolved_KTBins
  
} // end of namespace RESUM

#endif


