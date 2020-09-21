#ifndef RESUM_KT2_EE_H
#define RESUM_KT2_EE_H
#include "Analysis/ChannelAlgorithms/ChannelAlgorithm_Base.H"
#include <string>
#include <vector>
#include "ATOOLS/Math/Vector.H"


namespace RESUM {
  class KT2_ee : public ChannelAlgorithm_Base {
  public:
    KT2_ee(const ChAlg_Key& parameters);
    std::string Channel(const std::vector<ATOOLS::Vec4D>& ip,
                        const std::vector<ATOOLS::Flavour>& fl,
                        const size_t &nin,
                        std::vector<ATOOLS::Vec4D>* pout,
                        std::vector<ATOOLS::Flavour>* fout) override;

    std::vector<std::string> ChannelNames() override {return m_channelNames;}

    enum MODE {
               ALL = 0,
               BLAND = 1 << 0,
               ZDECAY = 1 << 1,
               BLAND_Z = ALL|BLAND|ZDECAY
    };
    
  private:
    
    double KT2(const ATOOLS::Vec4D &p1, const ATOOLS::Vec4D &p2,
               const std::vector<int>& fl1, const std::vector<int>& fl2, 
               int mode, int num_flavd) const;
    bool flavd(const std::vector<int>& fls) const;
    
    MODE m_mode = MODE::ALL;
    int m_nborn = -1;
    std::vector<std::string> m_channelNames;
    
  }; // end of class KT2_ee
  
} // end of namespace RESUM

#endif


