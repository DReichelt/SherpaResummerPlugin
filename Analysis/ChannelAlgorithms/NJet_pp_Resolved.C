#include "Analysis/ChannelAlgorithms/NJet_pp_Resolved.H"
#include "Analysis/ChannelAlgorithms/KT2_pp_Ordered.H"
#include "Observables/Algorithms/FastjetAlg.H"
#include "Analysis/Observable_Base.H"

using namespace RESUM;

NJet_pp_Resolved::NJet_pp_Resolved(const ChAlg_Key& parameters)
  : ChannelAlgorithm_Base(parameters) {
  int nmin = RESUM::to_type<int>(parameters.KwArg("NMIN","1"));
  int nmax = RESUM::to_type<int>(parameters.KwArg("NMAX"));
  std::vector<std::string> params = {"-1",parameters.KwArg("MODE","ALL")};
  for(int n=nmin; n<=nmax; n++) {
    params[0] = std::to_string(n);
    m_resolvers.emplace_back(new KT2_pp_Ordered({"KT2_pp_Ordered",params}));
    for(const std::string& name: m_resolvers.back()->ChannelNames(true)) {
      m_channelNames.push_back(name);
    }
    // m_channelNames.insert(m_channelNames.end(),
    //                       m_resolvers.back()->ChannelNames().begin(),
    //                       m_resolvers.back()->ChannelNames().end());
  }
}


std::string NJet_pp_Resolved::Channel(const std::vector<ATOOLS::Vec4D>& ip,
                                      const std::vector<ATOOLS::Flavour>& fl,
                                      const size_t &nin) {
  std::vector<ATOOLS::Vec4D> p;
  std::vector<ATOOLS::Flavour> f;
  for(int i=0; i<ip.size(); i++) {
    if(fl[i].Strong()) {
      p.push_back(ip[i]);
      f.push_back(fl[i]);
    }
  }
  const FJmaxPTjet fj(p,f,nin,Observable_Key("ChAlg",m_params));
  const int mult = fj.pseudoJets().size();
  return m_resolvers.at(mult-1)->Channel(ip,fl,nin,false);
}

DECLARE_GETTER(NJet_pp_Resolved,"NJet_pp_Resolved",ChAlg,ChAlg_Key);
ChAlg *ATOOLS::Getter<ChAlg,ChAlg_Key,NJet_pp_Resolved>::
operator()(const Parameter_Type &args) const 
{ return new NJet_pp_Resolved(args); }
void ATOOLS::Getter<ChAlg,ChAlg_Key,NJet_pp_Resolved>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"NJet_pp_Resolved"; }


