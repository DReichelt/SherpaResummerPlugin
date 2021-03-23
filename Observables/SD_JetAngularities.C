#include "Analysis/Observable_Base.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_Limits.H"
#include "FFunction/FFunctions.H"
#include <vector>       
#include <algorithm>       

#include "fastjet/tools/Recluster.hh"
#include "fastjet/contrib/SoftDrop.hh"


using namespace ATOOLS;

namespace RESUM {

  class SD_JetAngularities: public Observable_Base {
  public:

    SD_JetAngularities(const Observable_Key &args): 
    Observable_Base(args) {}

    Obs_Params Parameters
      (const std::vector<ATOOLS::Vec4D>& p,
       const std::vector<ATOOLS::Flavour>& fl,
       const size_t& l) {
      // TODO: this is a dummy
      return Obs_Params(2.0,0.0,0.0,0.0);
    }

    std::function<double(double,double&)> FFunction(const std::vector<ATOOLS::Vec4D>& p,
                                                    const std::vector<ATOOLS::Flavour>& fl, 
                                                    const RESUM::Params& params) {
      // TODO: this is a dummy
      return FFUNCTION::Additive;
    }


    virtual std::set<size_t> ResumMult() {return {3};}
    virtual size_t ResumQCDorderBorn() {return 1;};
    virtual size_t ResumQCDorderLO() {return 2;}
    virtual size_t ResumQCDorderNLO() {return 3;}
    


    double Value(const std::vector<Vec4D>& ip,
                 const std::vector<Flavour>& fl,
		 const size_t &nin) {

      if(ip.size() <= nin) return 0;
      msg_Debugging()<<"Start jet angularity.\n";
      // some parameters
      // TODO: make these settable
      const double R_jet = 0.5;          
      const double pt_min_cut = 0.;        
      const double alpha = 2.;
      
      // do this in runcard
      //pT_jet_cut = 500.; //GeV cut on the leading jet
      fastjet::Strategy strategy = fastjet::Best;
      fastjet::RecombinationScheme    recombScheme = fastjet::E_scheme;
      fastjet::JetDefinition* jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm,  R_jet, recombScheme, strategy);

      msg_Debugging()<<"Done fastjet settings.\n";
      
      fastjet::Recluster wta_recluster = fastjet::Recluster( fastjet::JetDefinition(fastjet::antikt_algorithm, 999.0, fastjet::WTA_pt_scheme) );

      msg_Debugging()<<"Done recluster.\n";
      
      std::vector<fastjet::PseudoJet> fjs(ip.size()-nin);
      for(size_t i=0; i<fjs.size(); i++) {
        msg_Debugging()<<"Add "<<ip.at(i+nin)<<"\n";
        fjs[i] = {ip.at(i+nin)};
      }

      fastjet::ClusterSequence cs(fjs, *jetDef);
      std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt( cs.inclusive_jets(pt_min_cut) );

      if (jets.size() < 1) {
        THROW(fatal_error, "Something went wrong!");
      }
      else {
        msg_Debugging()<<"Jets ok.\n";
      }
      // Run SoftDrop  
      //
      fastjet::contrib::SoftDrop sd_groomer(2.0, 0.05, 0.5);
      
      fastjet::PseudoJet sd_jet = sd_groomer(jets[0]);
      
      std::vector<fastjet::PseudoJet> sd_constituents = sd_jet.constituents();
      
      // This one we use  to compute dR
      fastjet::PseudoJet wta_jet = wta_recluster(jets[0]);    
      
      double lambda = 0.;
      
      for (const auto& el : sd_constituents) {
        
        double z = el.pt() / sd_jet.pt();
        
        double Delta = dR(el, wta_jet);
        
        lambda +=  z * pow(Delta / R_jet, alpha);
      }
      delete jetDef;
      return lambda;
    }           
    
    double dR (const fastjet::PseudoJet &parton, const fastjet::PseudoJet &wta_jet) { 
      //
      double dR = 0.;
      
      double delta_y = parton.rapidity() - wta_jet.rapidity();
      
      double delta_phi = parton.phi() - wta_jet.phi();
      
      dR = sqrt( delta_y * delta_y + delta_phi * delta_phi );
      
      return dR;
    }
    
    
  };// end of class Y1_II

}// end of namespace RESUM

using namespace RESUM;


DECLARE_GETTER(SD_JetAngularities,"SD_JetAngularities",Observable_Base,Observable_Key);
Observable_Base *ATOOLS::Getter<Observable_Base,Observable_Key,SD_JetAngularities>::
operator()(const Parameter_Type &args) const 
{ return new SD_JetAngularities(args); }
void ATOOLS::Getter<Observable_Base,Observable_Key,SD_JetAngularities>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"SD_JetAngularities"; }
