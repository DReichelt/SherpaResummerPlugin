#include "Observables/Algorithms/FastjetAlg.H"
#include "Analysis/Observable_Base.H"

using namespace RESUM;

FJmaxPTjet::FJmaxPTjet(const std::vector<ATOOLS::Vec4D>& p,
                       const std::vector<ATOOLS::Flavour>& fl,
                       const size_t &nin,
                       const Observable_Key& key)   {
  
  double R = to_type<double>(key.KwArg("R","0.5"));
  m_minPT = to_type<double>(key.KwArg("minPT","0"));
  double beta = to_type<double>(key.KwArg("sd_beta","2"));
  double zcut = to_type<double>(key.KwArg("sd_zcut","2"));
  fastjet::JetDefinition jetDef = fastjet::JetDefinition(fastjet::antikt_algorithm,  R, recombScheme, strategy);
  // TODO: get groom params from key
  fastjet::contrib::SoftDrop sd_groomer = {beta, zcut, R};          
  
  std::vector<fastjet::PseudoJet> fjs(p.size()-nin);
  for(size_t i=0; i<fjs.size(); i++) {
    if(fl[i+nin].Strong()) {
      msg_Debugging()<<"Add "<<p.at(i+nin)<<"\n";
      const ATOOLS::Vec4D& pi = p[i+nin];
      fjs[i] = {pi[1],pi[2],pi[3],pi[0]};
      fjs[i].set_user_index(i+nin);
    }
  }
  
  cs = {fjs, jetDef};
  
  if (pseudoJets().size() < 1) {
    THROW(fatal_error, "Something went wrong!");
  }
  else {
    msg_Debugging()<<"Jets ok.\n";
  }
  
  m_jet = pseudoJets()[0];
  m_sdJet = sd_groomer(m_jet);
  
  m_jetAxis = wta_recluster(m_jet);
  m_jets = {2,std::set<size_t>()};
  for(const fastjet::PseudoJet& j: m_jet.constituents() ) {
    m_jets[0].emplace(j.user_index());
  }
  for(const fastjet::PseudoJet& j: m_sdJet.constituents() ) {
    m_jets[1].emplace(j.user_index());
  }
  m_jetVectors = {{m_jet.E(),m_jet.px(),m_jet.py(),m_jet.pz()},
                  {m_sdJet.E(),m_sdJet.px(),m_sdJet.py(),m_sdJet.pz()}};
  m_jetAxes = {{m_jetAxis.px(),m_jetAxis.py(),m_jetAxis.pz()},
               {m_jetAxis.px(),m_jetAxis.py(),m_jetAxis.pz()}};
  m_jetScales = {m_jet.perp(), m_sdJet.perp()};
}
