#ifdef USING_FJCONTRIB 
#include "Analysis/ChannelAlgorithms/NJet_pp_Resolved_KTBins.H"
#include "Analysis/ChannelAlgorithms/KT2_pp_Ordered.H"
#include "Observables/Algorithms/FastjetAlg.H"
#include "Analysis/Observable_Base.H"
#include "Tools/StringTools.H"

using namespace RESUM;

NJet_pp_Resolved_KTBins::NJet_pp_Resolved_KTBins(const ChAlg_Key& parameters)
  : ChannelAlgorithm_Base(parameters) {
  int nmin = RESUM::to_type<int>(parameters.KwArg("NMIN","1"));
  int nmax = RESUM::to_type<int>(parameters.KwArg("NMAX"));
  m_edges = split(parameters.KwArg("EDGES"),"_");
  for(const std::string& e: m_edges) {
    m_binEdges.push_back(RESUM::to_type<double>(e));
  }
  std::vector<std::string> params = {"-1",parameters.KwArg("MODE","ALL"),"SUMNEUTRAL:"+parameters.KwArg("SUMNEUTRAL","0")};
  m_collapse = parameters.KwArg("COLLAPSE","0") != "0";
  for(int n=nmin; n<=nmax; n++) {
    params[0] = std::to_string(n);
    m_resolvers.emplace_back(new KT2_pp_Ordered({"KT2_pp_Ordered",params}));
    std::set<std::string> names;
    for(std::string name: m_resolvers.back()->ChannelNames(true)) {
      collapse(name);
      for(size_t i=0; i<m_edges.size()-1; i++) {
        names.insert(name+"_PtLead_"+m_edges[i]+"_"+m_edges[i+1]);
      }
      names.insert(name+"_PtLead_"+m_edges.back());
      collapse(name,"Other");
      for(size_t i=0; i<m_edges.size()-1; i++) {
        names.insert(name+"_PtLead_"+m_edges[i]+"_"+m_edges[i+1]);
      }
      names.insert(name+"_PtLead_"+m_edges.back());
      collapse(name,"EW");
      for(size_t i=0; i<m_edges.size()-1; i++) {
        names.insert(name+"_PtLead_"+m_edges[i]+"_"+m_edges[i+1]);
      }
      names.insert(name+"_PtLead_"+m_edges.back());
    }
    m_channelNames.assign(names.begin(), names.end());
    // m_channelNames.insert(m_channelNames.end(),
    //                       m_resolvers.back()->ChannelNames().begin(),
    //                       m_resolvers.back()->ChannelNames().end());
  }
}

std::string NJet_pp_Resolved_KTBins::collapse(std::string& name, std::string fl) {
  if(!m_collapse) return name;
  //if(name.find("Other") != std::string::npos) return name;
  if(fl=="") {
    fl = name.substr(name.find("TO")+2,1);
    if(fl=="g") fl = "Gluon";
    else if(fl=="q") fl = "Quark";
    else THROW(fatal_error,"Unknown flavour "+fl);
  }
  name = fl+name.substr(name.find("_"));
  return name;
}


std::string NJet_pp_Resolved_KTBins::Channel(const std::vector<ATOOLS::Vec4D>& ip,
                                             const std::vector<ATOOLS::Flavour>& fl,
                                             const size_t &nin,
                                             std::vector<ATOOLS::Vec4D>* pout,
                                             std::vector<ATOOLS::Flavour>* fout) {
  DEBUG_FUNC("");
  if(ip.size() < nin) {
    THROW(fatal_error,"Less particles than initial states.");
  }
  std::vector<ATOOLS::Vec4D> p = ip;
  std::vector<ATOOLS::Flavour> f = fl;
  // for(int i=0; i<ip.size(); i++) {
  //   if(fl[i].Strong()) {
  //     p.push_back(ip[i]);
  //     f.push_back(fl[i]);
  //     msg_Debugging()<<"Added "<<f.back()<<" "<<p.back()<<"\n";
  //   }
  // }
  const Observable_Key ChKey = {"ChAlg",m_params};
  int n = p.size()-nin;
  for(int i=0; i<ip.size(); i++) {
    if(not fl[i].Strong()) n--; 
  }
  if(n-1>m_resolvers.size()) 
    THROW(fatal_error, "No resolver for this multiplicity: "+std::to_string(n)+".");
  std::string channel = m_resolvers.at(n-1)->Channel(ip,fl,nin,false,&p,&f);
  // pt of leading (ungroomed) jet
  FJmaxPTjet fj(p,f,nin,ChKey);
  const double PT = fj.jetScales(0);
  bool ew = false;
  for(ATOOLS::Flavour flav: fj.apply(f,0)) {
    if(not flav.Strong()) {
      ew = true;
      collapse(channel,"EW");
      break;
    }
  }
  msg_Debugging()<<"Channel candidate: "<<channel<<".\n";
  if(not ew) {

    std::string flname = "";
    while(fj.pseudoJets().size() > 0 and fj.pseudoJets()[0].has_constituents() and fj.pseudoJets()[0].constituents().size() > 1) {
      for(auto& pd: fj.pseudoJets()[0].constituents()) {
        msg_Debugging()<<pd.e()<<" "<<pd.px()<<" "<<pd.py()<<" "<<pd.pz()<<"\n";
      }
      msg_Debugging()<<"Leading jet has "<<fj.pseudoJets()[0].constituents().size()<<" constituents, "<<p.size()-2<<" jets left to cluster.\n";
      n--;
      msg_Debugging()<<"Using cluster algorithm "<<n<<" of "<<m_resolvers.size()-1<<".\n";
      channel = m_resolvers.at(n-1)->Channel(ip,fl,nin,false,&p,&f);
      // create dummy flavours
      msg_Debugging()<<p<<" "<<f<<"\n";
      msg_Debugging()<<"New channel candidate: "<<channel<<"\n";
      fj = FJmaxPTjet(p,f,nin,ChKey);
    }
  }
  msg_Debugging()<<"Returning "<<channel<<".\n";
  if(pout) {
    *pout = p;
  }
  if(fout) {
    *fout = f;
  }
  if(PT > m_binEdges.back()) channel += "_PtLead_"+m_edges.back();
  else {
    for(size_t i=0; i<m_edges.size()-1; i++) {
      if(PT < m_binEdges[i+1] and PT > m_binEdges[i]) {
        channel += "_PtLead_"+m_edges[i]+"_"+m_edges[i+1];
      }
    }
  }
  if(m_collapse and not ew){
    std::string chname = "";
    if(fj.pseudoJets().size() > 0) {
      if(fj.apply(f,0)[0].IsQuark()) chname = "Quark";
      else chname = "Gluon";
    }
    else {
      chname = "Other";
    }
    collapse(channel,chname); 
  }
  return channel;//m_resolvers.at(mult-1)->Channel(ip,fl,nin,false);
}

DECLARE_GETTER(NJet_pp_Resolved_KTBins,"NJet_pp_Resolved_KTBins",ChAlg,ChAlg_Key);
ChAlg *ATOOLS::Getter<ChAlg,ChAlg_Key,NJet_pp_Resolved_KTBins>::
operator()(const Parameter_Type &args) const 
{ return new NJet_pp_Resolved_KTBins(args); }
void ATOOLS::Getter<ChAlg,ChAlg_Key,NJet_pp_Resolved_KTBins>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"NJet_pp_Resolved_KTBins"; }
#endif

