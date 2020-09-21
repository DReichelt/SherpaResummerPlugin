#ifndef RESUM_FJMAXPTJET_H
#define RESUM_FJMAXPTJET_H

#include "Algorithm.H"
#include "ATOOLS/Phys/Flavour.H"

#include "fastjet/tools/Recluster.hh"
#include "fastjet/contrib/SoftDrop.hh"


namespace RESUM {

  struct Observable_Key;
  
  class FJmaxPTjet : public Algorithm<double> {
  public:
    
    typedef std::shared_ptr<FJmaxPTjet> Ptr;
    
    FJmaxPTjet(const std::vector<ATOOLS::Vec4D>& p,
               const std::vector<ATOOLS::Flavour>& fl,
               const size_t &nin,
               const Observable_Key& key);

    std::vector<fastjet::PseudoJet> pseudoJets() const {return fastjet::sorted_by_pt( fastjet::SelectorAbsRapMax(m_maxRap)(cs.inclusive_jets(m_minPT)) );}

    const fastjet::PseudoJet LeadJet() const {return m_jet;}
    const fastjet::PseudoJet SDLeadJet() const {return m_sdJet;}
    const fastjet::PseudoJet LeadJetAxis() const {return m_jetAxis;}
    
  private:

    fastjet::ClusterSequence cs;
    fastjet::PseudoJet m_sdJet;
    fastjet::PseudoJet m_jet;
    fastjet::PseudoJet m_jetAxis;
    double m_minPT;
    double m_maxRap;
    // @TODO: make this dynamic
    fastjet::Strategy strategy = fastjet::Best;
    fastjet::RecombinationScheme recombScheme = fastjet::E_scheme;
    bool softDrop = true;
    fastjet::Recluster wta_recluster = fastjet::Recluster( fastjet::JetDefinition(fastjet::antikt_algorithm, 999.0, fastjet::WTA_pt_scheme) );

  };
}

#endif
