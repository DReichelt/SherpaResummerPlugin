#ifndef RESUM_KT2_PP_ORD_H
#define RESUM_KT2_PP_ORD_H
#include "Analysis/ChannelAlgorithms/ChannelAlgorithm_Base.H"
#include <string>
#include <vector>
#include "ATOOLS/Math/Vector.H"


namespace RESUM {
  class KT2_pp_Ordered : public ChannelAlgorithm_Base {
  public:
    KT2_pp_Ordered(const ChAlg_Key& parameters, bool orderByPT=true);
    std::string Channel(const std::vector<ATOOLS::Vec4D>& ip,
                        const std::vector<ATOOLS::Flavour>& fl,
                        const size_t &nin,
                        std::vector<ATOOLS::Vec4D>* pout,
                        std::vector<ATOOLS::Flavour>* fout) override;


    std::vector<std::string> ChannelNames() override {return m_channelNames;}

    enum MODE {
               ALL = 0,
               BLAND = 1      << 0,
               ZPROD = 1      << 1,
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
               const size_t n,
               const std::vector<ATOOLS::Vec4D> &neutralFinal) const;
    bool flavd(const std::vector<int>& fls) const;
    ATOOLS::Flavour afl(const std::vector<int>& fls) const {
      if(!flavd(fls)) return {21};
      ATOOLS::Flavour ret;
      for(int i=0; i<5; i++) {
        if(fls[i] != 0) {
          if(ret.IsDummy()) {
            ret = {fls[i]*(i+1)};
          }
          else {
            return ATOOLS::Flavour();
          }
        }
      }
      return ret;
    }

    MODE m_mode = MODE::ALL;
    int m_nborn = -1;
    bool m_sumNeutral = false;
    bool m_orderByPT = true;
    std::vector<std::string> m_channelNames;
    
  }; // end of class KT2_pp_Ordered
  
} // end of namespace RESUM

#endif


