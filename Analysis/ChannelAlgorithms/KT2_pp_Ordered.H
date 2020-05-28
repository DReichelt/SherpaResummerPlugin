#ifndef RESUM_KT2_PP_ORD_H
#define RESUM_KT2_PP_ORD_H
#include "Analysis/ChannelAlgorithms/ChannelAlgorithm_Base.H"
#include <string>
#include <vector>
#include "ATOOLS/Math/Vector.H"


namespace RESUM {
  class KT2_pp_Ordered : public ChannelAlgorithm_Base {
  public:
    KT2_pp_Ordered(const ChAlg_Key& parameters);
    std::string Channel(const std::vector<ATOOLS::Vec4D>& ip,
                        const std::vector<ATOOLS::Flavour>& fl,
                        const size_t &nin) override;

    std::vector<std::string> ChannelNames() override {return m_channelNames;}

    enum MODE {
               ALL = 0,
               BLAND = 1   << 0,
               ZPROD = 1   << 1,
               BLAND_Z = BLAND|ZPROD
    };
    
  private:
    
    double KT2(const ATOOLS::Vec4D &p1, const ATOOLS::Vec4D &p2,
               const std::vector<int>& fl1, const std::vector<int>& fl2, 
               int mode, int num_flavd) const;
    double KT2(const ATOOLS::Vec4D &p1,
               const double ktB2,
               const std::vector<int>& fl1,
               const std::vector<int>& fl2,
               int mode, int num_flavd) const;
    double KT2(const std::vector<ATOOLS::Vec4D> &ps,
               const double eta,
               const int beamIdx) const;
    double KT2(const std::vector<ATOOLS::Vec4D> &p,
               const double eta,
               const int beamIdx,
               const std::vector<int> imap,
               const int n) const;
    bool flavd(const std::vector<int>& fls) const;
    
    MODE m_mode = MODE::ALL;
    int m_nborn = -1;
    std::vector<std::string> m_channelNames;
    
  }; // end of class KT2_pp_Ordered
  
} // end of namespace RESUM

#endif


