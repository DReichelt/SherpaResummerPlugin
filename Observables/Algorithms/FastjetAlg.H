#ifndef RESUM_FJMAXPTJET_H
#define RESUM_FJMAXPTJET_H

#include "Algorithm.H"

#include "fastjet/tools/Recluster.hh"
#include "fastjet/contrib/SoftDrop.hh"


namespace RESUM {

  class FJmaxPTjet : public Algorithm<double> {
  public:
    
    typedef std::shared_ptr<FJmaxPTjet> Ptr;
    
    FJmaxPTjet(const std::vector<ATOOLS::Vec4D>& p,
               const std::vector<ATOOLS::Flavour>& fl,
               const size_t &nin)   {

      std::vector<fastjet::PseudoJet> fjs(p.size()-nin);
      for(size_t i=0; i<fjs.size(); i++) {
        msg_Debugging()<<"Add "<<p.at(i+nin)<<"\n";
        fjs[i] = {p.at(i+nin)};
        fjs[i].set_user_index(i);
      }

      cs = {fjs, jetDef};

      if (pseudoJets().size() < 1) {
        THROW(fatal_error, "Something went wrong!");
      }
      else {
        msg_Debugging()<<"Jets ok.\n";
      }

      m_jet = pseudoJets()[0];
      if(softDrop) {
        m_sdJet = sd_groomer(m_jet);
      }
      m_jetAxis = wta_recluster(m_jet);
      m_jets = {2,std::set<size_t>()};
      for(const fastjet::PseudoJet& j: m_jet.constituents() ) {
        m_jets[0].emplace(j.user_index());
      }
      for(const fastjet::PseudoJet& j: m_sdJet.constituents() ) {
        m_jets[1].emplace(j.user_index());
      }
      m_jetAxes = {{m_jetAxis.px(),m_jetAxis.py(),m_jetAxis.pz()},
                   {m_jetAxis.px(),m_jetAxis.py(),m_jetAxis.pz()}};
      m_Q2 = softDrop ? m_jet.pt() : m_sdJet.pt();
    }

    std::vector<fastjet::PseudoJet> pseudoJets() {return fastjet::sorted_by_pt( cs.inclusive_jets(pt_min) );}

    const fastjet::PseudoJet LeadJet() const {return m_jet;}
    const fastjet::PseudoJet SDLeadJet() const {return m_sdJet;}
    const fastjet::PseudoJet LeadJetAxis() const {return m_jetAxis;}
    
  private:

    fastjet::ClusterSequence cs;
    fastjet::PseudoJet m_sdJet;
    fastjet::PseudoJet m_jet;
    fastjet::PseudoJet m_jetAxis;
    
    // @TODO: make this dynamic
    fastjet::Strategy strategy = fastjet::Best;
    fastjet::RecombinationScheme recombScheme = fastjet::E_scheme;
    double R_jet = 0.5;
    double pt_min = 0;
    bool softDrop = true;
    fastjet::JetDefinition jetDef = fastjet::JetDefinition(fastjet::antikt_algorithm,  R_jet, recombScheme, strategy);
    fastjet::Recluster wta_recluster = fastjet::Recluster( fastjet::JetDefinition(fastjet::antikt_algorithm, 999.0, fastjet::WTA_pt_scheme) );
    fastjet::contrib::SoftDrop sd_groomer = {2.0, 0.05, 0.5};    
  };
}

#endif
