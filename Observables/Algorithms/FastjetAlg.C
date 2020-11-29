#ifdef USING_FJCONTRIB
#include "Observables/Algorithms/FastjetAlg.H"
#include "Analysis/Observable_Base.H"

using namespace RESUM;

FJmaxPTjet::FJmaxPTjet(const std::vector<ATOOLS::Vec4D>& p,
                       const std::vector<ATOOLS::Flavour>& fl,
                       const size_t &nin,
                       const Observable_Key& key)   {
  DEBUG_FUNC("");
  double R = to_type<double>(key.KwArg("R","0.5"));
  m_minPT = to_type<double>(key.KwArg("minPT","0"));
  m_maxRap = to_type<double>(key.KwArg("maxRap","-1"));
  if(m_maxRap < 0) m_maxRap = std::numeric_limits<double>::infinity();
  double maxAsym = to_type<double>(key.KwArg("maxAsym","-1"));
  if(maxAsym < 0) maxAsym = std::numeric_limits<double>::infinity();
  const double minDPhi = to_type<double>(key.KwArg("minDPhi","0"));
  const double minZPT = to_type<double>(key.KwArg("minZPT","0"));
  const double beta = to_type<double>(key.KwArg("beta","2"));
  const double zcut = to_type<double>(key.KwArg("zcut","0.05"));
  const bool ewCluster = std::set<std::string>({"no","NO","n","N","0"}).count(key.KwArg("ewCluster","no")) == 0;
  fastjet::JetDefinition jetDef = fastjet::JetDefinition(fastjet::antikt_algorithm,  R, recombScheme, strategy);
  // TODO: get groom params from key
  fastjet::contrib::SoftDrop sd_groomer = {beta, zcut, R};          


  std::vector<fastjet::PseudoJet> fjs;
  std::vector<fastjet::PseudoJet> strongfjs;
  std::vector<fastjet::PseudoJet> ewfjs;
  fastjet::PseudoJet Zfj = {0,0,0,0};
  for(size_t i=0; i<p.size()-nin; i++) {
    const ATOOLS::Vec4D& pi = p[i+nin];
    // msg_Out()<<"Add "<<fl[i+nin]<<" "<<p.at(i+nin)<<"\n";
    fjs.emplace_back(pi[1],pi[2],pi[3],pi[0]);
    fjs.back().set_user_index(i+nin);
    if(fl[i+nin].Strong()) {
      strongfjs.emplace_back(fjs.back());
      strongfjs.back().set_user_index(i+nin);
    }
    else {
      ewfjs.emplace_back(fjs.back());
      ewfjs.back().set_user_index(i+nin);
      Zfj += ewfjs.back();
    }
  }
  if (Zfj.pt() < minZPT) {
    msg_Debugging()<<"No viable Z.\n\n";
    m_jets.clear();
    m_jetVectors.clear();
    m_jetAxes.clear();
    m_jetScales.clear();
    return;
  }

  
  //cs = {fjs, jetDef};

  // m_pseudoJets = fastjet::sorted_by_pt((fastjet::SelectorAbsRapMax(m_maxRap)*fastjet::SelectorPtMin(m_minPT))(jetDef(strongfjs)));
  m_pseudoJets = fastjet::sorted_by_pt(fastjet::SelectorPtMin(m_minPT)(jetDef(strongfjs)));
  // m_pseudoJets.clear();
  // for(auto& fj: fastjet::sorted_by_pt( (fastjet::SelectorAbsRapMax(m_maxRap)*fastjet::SelectorPtMin(m_minPT))(jetDef(strongfjs)))) {
  //   bool add = true;
  //   //msg_Out()<<"Jet: ("<<fj.e()<<", "<<fj.px()<<", "<<fj.py()<<", "<<fj.pz()<<") {\n";
  //   for(auto& ffj: fj.constituents()) {
  //     //msg_Out()<<"\t\t"<<fl[ffj.user_index()]<<" "<<p[ffj.user_index()]<<"\n";
  //   }
  //   //msg_Out()<<"\t}\n";
  //   for(const auto& ewfj: ewfjs) {
  //     //msg_Out()<<ewfj.pt()<<" "<<fj.pt()<<" "<<fj.squared_distance(ewfj)<<" "<<ATOOLS::sqr(R)<<"\n";
  //     if(ewfj.pt() < 0.5*fj.pt()) //msg_Out()<<"Ignoring soft muon.\n";
  //     if(fj.squared_distance(ewfj) > ATOOLS::sqr(R)) //msg_Out()<<"No overlap.\n";
  //     add = add and (ewfj.pt() < 0.5*fj.pt() or fj.squared_distance(ewfj) > ATOOLS::sqr(R));
  //     add = add and fj.squared_distance(ewfj) > ATOOLS::sqr(R);
  //   }
  //   if(add) {
  //     //msg_Out()<<"Adding "<<"Jet: ("<<fj.e()<<", "<<fj.px()<<", "<<fj.py()<<", "<<fj.pz()<<")\n";
  //     m_pseudoJets.emplace_back(std::move(fj));
  //   }
  //   else {
  //     //msg_Out()<<"Vetoing "<<"Jet: ("<<fj.e()<<", "<<fj.px()<<", "<<fj.py()<<", "<<fj.pz()<<")\n";
  //   }
  // }
  // //msg_Out()<<"\n\n";
  

  if (pseudoJets().size() < 1) {
    msg_Debugging()<<"No viable jets.\n\n";
    m_jets.clear();
    m_jetVectors.clear();
    m_jetAxes.clear();
    m_jetScales.clear();
  }
  else {
    msg_Debugging()<<"Jets ok.\n\n";
    m_jet = pseudoJets()[0];
    if (std::abs(m_jet.rapidity()) > m_maxRap or
        std::abs(Zfj.delta_phi_to(m_jet)) < minDPhi or 
        std::abs((m_jet.pt() - Zfj.pt()) / (m_jet.pt()+Zfj.pt())) > maxAsym) {
      msg_Debugging()<<"Cut on Z or jet failed.\n\n";
      m_jets.clear();
      m_jetVectors.clear();
      m_jetAxes.clear();
      m_jetScales.clear();
      return;
    }

    // std::vector<fastjet::PseudoJet> altJets;
    // for(auto& fj: fastjet::sorted_by_pt( (fastjet::SelectorAbsRapMax(m_maxRap)*fastjet::SelectorPtMin(m_minPT))(jetDef(fjs))) ) {
    //   bool add = true;
    //   //msg_Out()<<"Jet: ("<<fj.e()<<", "<<fj.px()<<", "<<fj.py()<<", "<<fj.pz()<<") {\n";
    //   for(auto& ffj: fj.constituents()) {
    //     //msg_Out()<<"\t\t"<<fl[ffj.user_index()]<<" "<<p[ffj.user_index()]<<"\n";
    //   }
    //   //msg_Out()<<"\t}\n";
    //   for(const auto& ewfj: ewfjs) {
    //     //msg_Out()<<ewfj.pt()<<" "<<fj.pt()<<" "<<fj.squared_distance(ewfj)<<" "<<ATOOLS::sqr(R)<<"\n";
    //     if(ewfj.pt() < 0.5*fj.pt()) //msg_Out()<<"Ignoring soft muon.\n";
    //     if(fj.squared_distance(ewfj) > ATOOLS::sqr(R)) //msg_Out()<<"No overlap.\n";
    //     add = add and (ewfj.pt() < 0.5*fj.pt() or fj.squared_distance(ewfj) > ATOOLS::sqr(R));
    //     add = add and fj.squared_distance(ewfj) > ATOOLS::sqr(R);
    //   }
    //   if(add) {
    //     //msg_Out()<<"Adding "<<"Jet: ("<<fj.e()<<", "<<fj.px()<<", "<<fj.py()<<", "<<fj.pz()<<")\n";
    //     altJets.emplace_back(std::move(fj));
    //   }
    //   else {
    //       //msg_Out()<<"Vetoing "<<"Jet: ("<<fj.e()<<", "<<fj.px()<<", "<<fj.py()<<", "<<fj.pz()<<")\n";
    //     }
    // }
    // //msg_Out()<<"\n\n";
    // ////msg_Out()<<altJets[0].constituents().size()<<" "<<m_jet.constituents().size()<<"\n";
    // if(altJets.size() < 1) {
    //     //msg_Out()<<"Jet: ("<<m_jet.e()<<", "<<m_jet.px()<<", "<<m_jet.py()<<", "<<m_jet.pz()<<") {\n";
    //     for(auto& ffj: m_jet.constituents()) {
    //       //msg_Out()<<"\t\t"<<fl[ffj.user_index()]<<" "<<p[ffj.user_index()]<<"\n";
    //     }
    //     msg_Out()<<"No jet from weak clustering.\n";
    //     return;
    //     THROW(fatal_error,"");
    // }
    // if(altJets[0].constituents().size() != m_jet.constituents().size()) {
    //   msg_Out()<<"Jet content changed due to ew clustering.\n";
    //   msg_Out()<<"Jet: ("<<m_jet.e()<<", "<<m_jet.px()<<", "<<m_jet.py()<<", "<<m_jet.pz()<<") {\n";
    //     for(auto& ffj: m_jet.constituents()) {
    //       msg_Out()<<"\t\t"<<fl[ffj.user_index()]<<" "<<p[ffj.user_index()]<<"\n";
    //     }
    //     msg_Out()<<"Alt. Jet: ("<<altJets[0].e()<<", "<<altJets[0].px()<<", "<<altJets[0].py()<<", "<<altJets[0].pz()<<") {\n";
    //     for(auto& ffj: altJets[0].constituents()) {
    //       msg_Out()<<"\t\t"<<fl[ffj.user_index()]<<" "<<p[ffj.user_index()]<<"\n";
    //     }
    //     return;
    //     THROW(fatal_error,"");
    // }
    
    m_sdJet = sd_groomer(m_jet);

    // for(auto& fj: m_jet.constituents()) {

    //   if(fl[fj.user_index()].Strong() and m_jet.constituents().size()>1) {
    //     msg_Out()<<"EW Jet: ("<<m_jet.e()<<", "<<m_jet.px()<<", "<<m_jet.py()<<", "<<m_jet.pz()<<") {\n";
    //     for(auto& ffj: m_jet.constituents()) {
    //       msg_Out()<<"\t\t"<<fl[ffj.user_index()]<<" "<<p[ffj.user_index()]<<"\n";
    //     }
    //     msg_Out()<<"\t}\n";
    //     msg_Out()<<"Distance to Z: "<<Zfj.delta_phi_to(pseudoJets()[0])<<"\n";
    //     msg_Out()<<"Z: ("<<m_jet.e()<<", "<<m_jet.px()<<", "<<m_jet.py()<<", "<<m_jet.pz()<<") {\n";
    //     for(auto& ffj: ewfjs) {
    //       msg_Out()<<"\t\t"<<fl[ffj.user_index()]<<" "<<p[ffj.user_index()]<<"\n";
    //     }
    //     msg_Out()<<"\t}\n";
    //     msg_Out()<<"\n\n";
    //     // THROW(fatal_error,"");
    //   }
    // }
    

    m_jetAxis = wta_recluster(m_jet);
    m_jets = {2,std::set<size_t>()};
    if(m_jet.has_constituents()) {
      for(const fastjet::PseudoJet& j: m_jet.constituents() ) {
        m_jets[0].emplace(j.user_index());
      }
    }
    else {
      m_jets[0].emplace(m_jet.user_index());
    }
    if(m_sdJet.has_constituents()) {
      for(const fastjet::PseudoJet& j: m_sdJet.constituents() ) {
        m_jets[1].emplace(j.user_index());
      }
    }
    else {
      m_jets[1].emplace(m_sdJet.user_index());
    }
    m_jetVectors = {{m_jet.E(),m_jet.px(),m_jet.py(),m_jet.pz()},
                    {m_sdJet.E(),m_sdJet.px(),m_sdJet.py(),m_sdJet.pz()}};
    m_jetAxes = {{m_jetAxis.px(),m_jetAxis.py(),m_jetAxis.pz()},
                 {m_jetAxis.px(),m_jetAxis.py(),m_jetAxis.pz()}};
    m_jetScales = {m_jet.perp(), m_sdJet.perp()};
  }  
}
#endif
