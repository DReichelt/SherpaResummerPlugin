#include "Analysis/ChannelAlgorithms/NJet_pp_Resolved.H"
#include "Analysis/ChannelAlgorithms/KT2_pp_Ordered.H"
#include "Observables/Algorithms/FastjetAlg.H"
#include "Analysis/Observable_Base.H"

using namespace RESUM;

NJet_pp_Resolved::NJet_pp_Resolved(const ChAlg_Key& parameters)
  : ChannelAlgorithm_Base(parameters) {
  int nmin = RESUM::to_type<int>(parameters.KwArg("NMIN","1"));
  int nmax = RESUM::to_type<int>(parameters.KwArg("NMAX"));
  std::vector<std::string> params = {"-1",parameters.KwArg("MODE","ALL"),"SUMNEUTRAL:"+parameters.KwArg("SUMNEUTRAL","0")};
  m_collapse = parameters.KwArg("COLLAPSE","0") != "0";
  for(int n=nmin; n<=nmax; n++) {
    params[0] = std::to_string(n);
    m_resolvers.emplace_back(new KT2_pp_Ordered({"KT2_pp_Ordered",params}));
    std::set<std::string> names;
    for(std::string name: m_resolvers.back()->ChannelNames(true)) {
      names.insert(collapse(name));
      //m_channelNames.push_back(collapse(name));
    }
    m_channelNames.assign(names.begin(), names.end());
    // m_channelNames.insert(m_channelNames.end(),
    //                       m_resolvers.back()->ChannelNames().begin(),
    //                       m_resolvers.back()->ChannelNames().end());
  }
}

std::string NJet_pp_Resolved::collapse(std::string& name) {
  if(!m_collapse) return name;
  msg_Out()<<name<<" -> ";
  std::string fl = name.substr(name.find("TO")+3,name.find("TO")+4);
  if(fl=="g") fl = "Gluon";
  else if(fl=="q") fl = "Quark";
  else THROW(fatal_error,"Unknown flavour "+fl);
  name = fl+name.substr(name.find("_"));
  msg_Out()<<name<<"\n";
  return name;
}


std::string NJet_pp_Resolved::Channel(const std::vector<ATOOLS::Vec4D>& ip,
                                      const std::vector<ATOOLS::Flavour>& fl,
                                      const size_t &nin,
                                      std::vector<ATOOLS::Vec4D>* pout) {
  DEBUG_FUNC("");
  std::vector<ATOOLS::Vec4D> p;
  std::vector<ATOOLS::Flavour> f;
  for(int i=0; i<ip.size(); i++) {
    if(fl[i].Strong()) {
      p.push_back(ip[i]);
      f.push_back(fl[i]);
      msg_Debugging()<<"Added "<<f.back()<<" "<<p.back()<<"\n";
    }
  }
  
  size_t n = p.size()-2;
  if(n-1>m_resolvers.size()) 
    THROW(fatal_error, "No resolver for this multiplicity: "+std::to_string(n)+".");
  std::string channel = m_resolvers.at(n-1)->Channel(ip,fl,nin,false,&p);
  msg_Debugging()<<"Channel candidate: "<<channel<<".\n";
  FJmaxPTjet fj(p,f,nin,Observable_Key("ChAlg",m_params));
  while(fj.pseudoJets()[0].has_constituents() and fj.pseudoJets()[0].constituents().size() > 1) {
    for(auto& p: fj.pseudoJets()[0].constituents()) {
      msg_Debugging()<<p.e()<<" "<<p.px()<<" "<<p.py()<<" "<<p.pz()<<"\n";
    }
    msg_Debugging()<<"Leading jet has "<<fj.pseudoJets()[0].constituents().size()<<" constituents, "<<p.size()-2<<" jets left to cluster.\n";
    n--;
    msg_Debugging()<<"Using cluster algorithm "<<n<<" of "<<m_resolvers.size()-1<<".\n";
    channel = m_resolvers.at(n-1)->Channel(ip,fl,nin,false,&p);
    f = std::vector<ATOOLS::Flavour>(p.size(),{21});
    msg_Debugging()<<p<<" "<<f<<"\n";
    msg_Debugging()<<"New channel candidate: "<<channel<<"\n";
    fj = FJmaxPTjet(p,f,nin,Observable_Key("ChAlg",m_params));
  }
  // std::vector<ATOOLS::Flavour> lead = fj.apply(f,0);
  
  // if(f.size()==5 and mult<3) {
  //   for(int i=0; i<ip.size(); i++) {
  //     if(fl[i].Strong()) {
  //       msg_Out()<<ip[i]<<" "<<fl[i]<<"\n";
  //     }
  //   }
  //   msg_Out()<<"Found "<<mult<<" jets.\n";
  //   for(auto& f: lead) msg_Out()<<f<<" ";
  //   msg_Out()<<" -> "<<channel<<"\n\n\n";
  // }
  collapse(channel);
  msg_Debugging()<<"Returning "<<channel<<".\n";
  if(pout) {
    *pout = p;
  }
  return channel;//m_resolvers.at(mult-1)->Channel(ip,fl,nin,false);
}

DECLARE_GETTER(NJet_pp_Resolved,"NJet_pp_Resolved",ChAlg,ChAlg_Key);
ChAlg *ATOOLS::Getter<ChAlg,ChAlg_Key,NJet_pp_Resolved>::
operator()(const Parameter_Type &args) const 
{ return new NJet_pp_Resolved(args); }
void ATOOLS::Getter<ChAlg,ChAlg_Key,NJet_pp_Resolved>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"NJet_pp_Resolved"; }


