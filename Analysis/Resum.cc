// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalPartons.hh"
#include "Rivet/Projections/Beam.hh"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "Main/Resum.H"
#include "Tools/StringTools.H"
#include "Observable_Base.H"

namespace Rivet {


  /// @brief Add a short analysis description here
  class Resum : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(Resum);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      ATOOLS::Data_Reader reader(",",";","!","=");
      ATOOLS::Algebra_Interpreter* ip = reader.Interpreter();
      p_resum=(RESUM::Resum*)ATOOLS::ToType<void*>(ATOOLS::rpa->gen.Variable("SHOWER_GENERATOR"));
      if (dynamic_cast<RESUM::Resum*>(p_resum) == nullptr) {
        THROW(fatal_error,"Resummer plugin not loaded");
      }
      p_resum->ResetObservables();
      string params = reader.GetValue<string>("RESUM::Observables","NO");
      if(params=="NO") THROW(fatal_error,"No Observables defined");
      size_t begin = params.find("Obs[");
      if(begin != string::npos) begin += 4;
      size_t end = params.find("]")-begin;
      while(begin != string::npos) {
        //cout<<"New Obs: ";
        const string& obs = params.substr(begin,end);
        //cout<<obs<<endl;
        params = params.substr(begin+end+1,params.size()-begin-end-1);
        //cout<<"Remaining: "<<params<<endl;
        begin = params.find("Obs[");
        if(begin != string::npos) begin += 4;
        end = params.find("]")-begin;
        nObs++;
        vector<string> ps = RESUM::split(obs,"&");
        string obsName = ps[0];
        ObsNames.push_back(obsName);
        double xmin = ATOOLS::ToType<double>(ip->Interprete(ps[1]));
        double xmax = ATOOLS::ToType<double>(ip->Interprete(ps[2]));
        double nbin = ATOOLS::ToType<double>(ip->Interprete(ps[3]));
        string htype = ps[4];
        //cout<<obsName<<" "<<xmin<<" "<<xmax<<" "<<nbin<<"\n";
        _h_NLL.push_back(Histo1D(nbin,xmin,xmax));
        _h_expLO.push_back(Histo1D(nbin,xmin,xmax));
        _h_expLO_diff.push_back(bookHisto1D(obsName+"expLO_diff",nbin,xmin,xmax));
        _h_expNLO.push_back(Histo1D(nbin,xmin,xmax));
        _s_NLL.push_back(bookScatter2D(obsName+"_NLL"));
        _s_expLO.push_back(bookScatter2D(obsName+"_expLO"));
        _s_expNLO.push_back(bookScatter2D(obsName+"_expNLO"));
        vector<double> xvals(_h_NLL.back().xEdges().size());
        if(htype == "LinErr")
          for(size_t i=0; i<xvals.size(); i++)
            xvals[i] = _h_NLL.back().xEdges()[i];
        if(htype == "LogErr")
          for(size_t i=0; i<xvals.size(); i++)
            xvals[i] = pow(10,_h_NLL.back().xEdges()[i]);
        if(htype == "LnErr")
          for(size_t i=0; i<xvals.size(); i++)
            xvals[i] = exp(_h_NLL.back().xEdges()[i]);
        p_resum->AddObservable(obsName,xvals);
      }
      declare(FinalPartons(), "PARTONS");
      declare(Beam(), "BEAM");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // RESUM::Observable_Base::Ptr tst = nullptr;
      // tst.reset(RESUM::Observable_Getter::GetObject("Y3_Durham",RESUM::Observable_Key("Y3_Durham")));
      // const Particles& ps = apply<FinalPartons>(event, "PARTONS").particles();
      // const ParticlePair& in = apply<Beam>(event, "BEAM").beams();
      // size_t mult = ps.size()+2;
      // vector<ATOOLS::Vec4D> moms(mult);
      // vector<ATOOLS::Flavour> flavs(mult);
      // moms[0] = momconv(in.first);
      // moms[1] = momconv(in.second);
      // flavs[0] = in.first.pdgId();
      // flavs[1] = in.second.pdgId();
      // for(size_t i=2; i<mult; i++) {
      //   moms[i] = momconv(ps[i-2]);
      //   flavs[i] = ps[i-2].pdgId();
      // }
      // if(tst->Value(moms,flavs) < 0.008) cout<<"Y3 = "<<tst->Value(moms,flavs)<<endl;
      
      if(!p_resum) THROW(fatal_error, "Plugin not loaded.");
      if(!p_resum->Initialized()) return;
      const double weight = event.weight();
      const double ntrial = event.genEvent()->weights()[3];
      N += ntrial;
      sumW += weight;      
      string channel = "";
      for(Particle p: apply<FinalPartons>(event,"PARTONS").particles()) {
        // @TODO: make this work if other particles are present!
        if(p.pid() == 21) channel = "g"+channel;
        else channel += "q";
      }
      if(_h_NLL_channels.find(channel) == _h_NLL_channels.end()) {
        _h_NLL_channels[channel] = vector<Histo1D>(_h_NLL.size());
        _s_NLL_channels[channel] = vector<Scatter2DPtr>(_h_NLL.size());
        _h_expLO_channels[channel] = vector<Histo1D>(_h_expLO.size());
        _s_expLO_channels[channel] = vector<Scatter2DPtr>(_h_expLO.size());
        _h_expNLO_channels[channel] = vector<Histo1D>(_h_expNLO.size());
        _s_expNLO_channels[channel] = vector<Scatter2DPtr>(_h_expNLO.size());
        for(size_t i=0; i<_h_NLL.size(); i++) {
          _h_NLL_channels[channel][i] = Histo1D(_h_NLL[i]);
          _h_NLL_channels[channel][i].reset();
          _s_NLL_channels[channel][i] = bookScatter2D(ObsNames[i]+"_NLL_Channel_"+channel);
          _h_expLO_channels[channel][i] = Histo1D(_h_expLO[i]);
          _h_expLO_channels[channel][i].reset();
          _s_expLO_channels[channel][i] = bookScatter2D(ObsNames[i]+"_expLO_Channel_"+channel);
          _h_expNLO_channels[channel][i] = Histo1D(_h_expNLO[i]);
          _h_expNLO_channels[channel][i].reset();
          _s_expNLO_channels[channel][i] = bookScatter2D(ObsNames[i]+"_expNLO_Channel_"+channel);
        }
      }
      for(size_t k=0; k<nObs; k++) {
        const double x = _h_NLL[k].xMin()-1.;
        _h_NLL[k].fill(x,p_resum->m_resNLL[k][0]*weight);
        _h_NLL_channels[channel][k].fill(x,p_resum->m_resNLL[k][0]*weight);
        _h_expLO[k].fill(x,p_resum->m_resExpLO[k][0]*weight);
        _h_expLO_channels[channel][k].fill(x,p_resum->m_resExpLO[k][0]*weight);
        _h_expNLO[k].fill(x,p_resum->m_resExpNLO[k][0]*weight);
        _h_expNLO_channels[channel][k].fill(x,p_resum->m_resExpNLO[k][0]*weight);
        for(size_t i=0; i<_h_NLL[k].bins().size(); i++) {
          const double x = _h_NLL[k].bins()[i].xMid();
          _h_NLL[k].fill(x,p_resum->m_resNLL[k][i+1]*weight);
          _h_NLL_channels[channel][k].fill(x,p_resum->m_resNLL[k][i+1]*weight);
          _h_expLO[k].fill(x,p_resum->m_resExpLO[k][i+1]*weight);
          _h_expLO_channels[channel][k].fill(x,p_resum->m_resExpLO[k][i+1]*weight);
          _h_expNLO[k].fill(x,p_resum->m_resExpNLO[k][i+1]*weight);
          _h_expNLO_channels[channel][k].fill(x,p_resum->m_resExpNLO[k][i+1]*weight);
        }
        
        // const std::pair<int, double>& res = p_resum->Result(k);
        for(size_t i=0; i<_h_expLO_diff[k]->bins().size()-1; i++) {

          // if(i==res.first) {
          //   // cout<<"Results: "<<res.second<<" "<<(p_resum->m_resExpLO[k][i+1]-p_resum->m_resExpLO[k][i])<<endl;
          //   const double x = _h_expLO_diff[k]->bins()[i].xMid();
          // _h_expLO_diff[k]->fill(x,(p_resum->m_resExpLO[k][i+1]-p_resum->m_resExpLO[k][i])*weight);

          // }
          // // else {
          // //   cout<<(p_resum->m_resExpLO[k][i+1]-p_resum->m_resExpLO[k][i])<<endl;
          // // }
          const double x = _h_expLO_diff[k]->bins()[i].xMid();
          _h_expLO_diff[k]->fill(x,(p_resum->m_resExpLO[k][i+1]-p_resum->m_resExpLO[k][i])*weight);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      cout<<sumOfWeights()/numEvents()<<endl;
      cout<<sumW/N<<endl;
      cout<<crossSection()<<endl;
      cout<<sumOfWeights()<<" "<<numEvents()<<endl;
      // const double xs = sumOfWeights()/numEvents();
      for(size_t k=0; k<nObs; k++) {
        //const double xs = _h_NLL[k].bins().back().sumW()/(numEvents()-1.);
        const int evts = N;//round(xs/crossSection()*(numEvents()-1.));
        cout<<evts<<" "<<N<<" "<<numEvents()<<endl;
        

        _h_expLO_diff[k]->scaleW(1./evts);
        
        _h_NLL[k].scaleW(1./(evts));
        _s_NLL[k]->addPoint(_h_NLL[k].xMin(),_h_NLL[k].underflow().sumW());
        for(auto& b: _h_NLL[k].bins())
          _s_NLL[k]->addPoint(b.xMax(),b.sumW());
        for(auto& ch: _h_NLL_channels) {
          ch.second[k].scaleW(1./evts);
          _s_NLL_channels[ch.first][k]->addPoint(ch.second[k].xMin(),ch.second[k].underflow().sumW());
          for(auto& b: ch.second[k].bins())
            _s_NLL_channels[ch.first][k]->addPoint(b.xMax(),b.sumW());
        }

        _h_expLO[k].scaleW(1./(evts));
        _s_expLO[k]->addPoint(_h_expLO[k].xMin(),_h_expLO[k].underflow().sumW());
        for(auto& b: _h_expLO[k].bins())
          _s_expLO[k]->addPoint(b.xMax(),b.sumW());
        for(auto& ch: _h_expLO_channels) {
          ch.second[k].scaleW(1./evts);
          _s_expLO_channels[ch.first][k]->addPoint(ch.second[k].xMin(),ch.second[k].underflow().sumW());
          for(auto& b: ch.second[k].bins())
            _s_expLO_channels[ch.first][k]->addPoint(b.xMax(),b.sumW());
        }

        _h_expNLO[k].scaleW(1./(evts));
        _s_expNLO[k]->addPoint(_h_expNLO[k].xMin(),_h_expNLO[k].underflow().sumW());
        for(auto& b: _h_expNLO[k].bins()) {
          _s_expNLO[k]->addPoint(b.xMax(),b.sumW());
        }
        for(auto& ch: _h_expNLO_channels) {
          ch.second[k].scaleW(1./evts);
          _s_expNLO_channels[ch.first][k]->addPoint(ch.second[k].xMin(),ch.second[k].underflow().sumW());
          for(auto& b: ch.second[k].bins()) {
            _s_expNLO_channels[ch.first][k]->addPoint(b.xMax(),b.sumW());
          }
        }

      }
    }

    //@}


    /// @name Histograms
    //@{
    size_t nObs = 0;
    vector<string> ObsNames;
    double N = 0;
    double sumW = 0;
    vector<Histo1D> _h_NLL;
    vector<Histo1D> _h_expLO;
    vector<Histo1DPtr> _h_expLO_diff;
    vector<Histo1D> _h_expNLO;
    map<string,vector<Histo1D>> _h_NLL_channels;
    map<string,vector<Histo1D>> _h_expLO_channels;
    map<string,vector<Histo1D>> _h_expNLO_channels;
    vector<Scatter2DPtr> _s_NLL;
    vector<Scatter2DPtr> _s_expLO;
    vector<Scatter2DPtr> _s_expNLO;
    map<string,vector<Scatter2DPtr>> _s_NLL_channels;
    map<string,vector<Scatter2DPtr>> _s_expLO_channels;
    map<string,vector<Scatter2DPtr>> _s_expNLO_channels;

    //@}

  private:
    RESUM::Resum *p_resum;

  private:
    ATOOLS::Vec4D momconv(const FourMomentum& other) {
      return {other.E(),other.px(),other.py(),other.pz()};
    }

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(Resum);


}
