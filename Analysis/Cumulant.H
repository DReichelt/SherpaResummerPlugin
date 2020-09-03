#ifndef RESUM_CUMULANT_H
#define RESUM_CUMULANT_H
#include <vector>
#include <valarray>
#include <memory>
#include <numeric>

#include "ATOOLS/Org/Gzip_Stream.H"


namespace RESUM {

  class Cumulant {
  
  public:

    typedef std::shared_ptr<Cumulant> Ptr;
    typedef std::map<std::string,std::valarray<size_t>> VarWeightMap;

    static std::valarray<double> Edges(const std::string& type, double lower, 
                                       double upper, int nbin) {
      if (type=="LogErr") return Edges(11,lower,upper,nbin);
      if (type=="LnErr") return Edges(21,lower,upper,nbin);
      if (type=="LinErr") return Edges(1,lower,upper,nbin);
    }

    static std::valarray<double> Edges(int type, double lower, 
                                       double upper, int nbin) {
      std::valarray<double> ret(nbin+1);
      switch(type) {
      case 1:
        ret[0] = lower;
        ret[nbin] = upper;
        break;
      case 11:
        ret[0] = std::log10(lower);
        ret[nbin-1] = std::log10(upper);
        break;
      case 21:
        ret[0] = std::log(lower);
        ret[nbin] = std::log(upper);
        break;
      default:
        THROW(fatal_error, "Histo type unknown.");
      }
      const double binsize = (ret[nbin-1]-ret[0])/nbin;
      for(int i=1; i<nbin; i++) ret[i] = ret[i-1]+binsize;
      switch(type) {
      case 1:
        break;
      case 11:
        ret = std::pow(10.,ret);
        break;
      case 21:
        ret = std::exp(ret);
        break;
      }
      return ret;
    }
    
    // Cumulant(std::valarray<double> edges, const std::string& name, int oform=1) 
    //   : m_name(name), m_edges(edges), m_nEdges(m_edges.size()), m_oform(oform),  
    //     m_weights(0.,m_nEdges), m_weight2s(0.,m_nEdges),
    //     m_barWeights(0.,m_nEdges), m_barWeight2s(0.,m_nEdges) {}

    Cumulant(std::valarray<double> edges, const std::string& name, 
             std::shared_ptr<VarWeightMap> varweights=nullptr, int oform=1)
      : m_name(name), m_edges(edges), m_nEdges(m_edges.size()), m_oform(oform),  
        p_varweights(varweights),
        m_nVarWeights((p_varweights?p_varweights->size():1)),
        m_nWeights(m_nEdges*m_nVarWeights),
        m_weights(0.,m_nWeights),
        m_weight2s(0.,m_nWeights),
        m_barWeights(0.,m_nWeights), 
        m_barWeight2s(0.,m_nWeights)  
    {
      DEBUG_FUNC(name);
      // if(!p_varweights) {
      //   msg_Out()<<"Initialized Cumulant without varweight.\n";
      //   if(varweights)
      //     msg_Out()<<"But had some earlier.\n";
      // }
      // else{
      //   msg_Out()<<"Initialize Cumulant with varweight.\n";
      // }
    }


    Cumulant& operator+=(Cumulant other) {
      if(m_nEdges == other.m_nEdges) {
        m_weights += other.m_weights;
        m_weight2s += other.m_weight2s;
        m_barWeights += other.m_barWeights;
        m_barWeight2s += other.m_barWeight2s;
        m_fills += other.m_fills;
      }
      else {
        THROW(fatal_error,"Merging of cumulants not implemented.");
      }
      return *this;
    }

    void Fill(double val, double weight, int ncount) {
      DEBUG_FUNC(val);
      Fill(val,weight,nullptr,ncount);
    }

    void FillMCB(double val, double weight, int ncount) {
      DEBUG_FUNC(val);
      FillMCB(val,weight,nullptr,ncount);
    }

    void Fill(double val, double cweight, 
              std::map<std::string,double>* weights, int ncount) {
      DEBUG_FUNC(val);
      if(!p_varweights) {
        InitVarWeights(weights);
        DO_STACK_TRACE;
        exit(1);
      }
      if(!weights) Fill(val,{{"CENTRAL",cweight}},ncount);
      else {
        (*weights)["CENTRAL"] = cweight;
        Fill(val,weights,ncount);
      }
    }

    void FillMCB(double val, double cweight, 
                 std::map<std::string,double>* weights, int ncount) {
      DEBUG_FUNC(val);
      if(!p_varweights) {
        InitVarWeights(weights);
        DO_STACK_TRACE;
        exit(1);
      }
      if(!weights) FillMCB(val,{{"CENTRAL",cweight}},ncount);
      else {
        (*weights)["CENTRAL"] = cweight;
        FillMCB(val,weights,ncount);
      }      
    }


    void Fill(double val, std::map<std::string,double> weights, int ncount) {
      DEBUG_FUNC(val);
      Fill(val,&weights,ncount);
    }

    void FillMCB(double val, std::map<std::string,double> weights, int ncount) {
      DEBUG_FUNC(val);
      FillMCB(val,&weights,ncount);
    }



    void Fill(double val, std::map<std::string,double>* weights, int ncount) {
      DEBUG_FUNC(val);
      m_fills += ncount;
      if(!p_varweights) {
        InitVarWeights(weights);
        DO_STACK_TRACE;
        exit(1);
      }
      if(!weights) return;
      for(const auto& var: *p_varweights) {
        const double w = (*weights)[var.first];
        // for(auto f: m_edges) msg_Out()<<f<<" ";
        // msg_Out()<<"\n";
        const std::valarray<bool> fill = m_edges > val;
        const std::valarray<bool> barFill = not fill;
        // for(auto f: fill) msg_Out()<<f<<" ";
        // msg_Out()<<"\n";
        // for(auto f: var.second[fill]) msg_Out()<<f<<" ";
        // msg_Out()<<"\n";
        const size_t nfills = std::count(std::begin(fill),std::end(fill),true);
        std::valarray<double> wBarFill = std::valarray<double>(w,m_nEdges-nfills); 
        std::valarray<double> wFill = std::valarray<double>(w,nfills); 
        if(ncount > 0) {
          m_weights[var.second[fill]] += wFill;
          m_barWeights[var.second[barFill]] += wBarFill;
          wFill *= wFill;
          wBarFill *= wBarFill;
          m_weight2s[var.second[fill]] += std::move(wFill);
          m_barWeight2s[var.second[barFill]] += std::move(wBarFill);
        }
        else {
          m_weights[var.second[fill]] += std::move(wFill);
          m_barWeights[var.second[barFill]] += std::move(wBarFill);
        }
      }
    }

    void FillMCB(double val, std::map<std::string,double>* weights, int ncount) {
      DEBUG_FUNC(val);
      if(!p_varweights) {
        InitVarWeights(weights);
        DO_STACK_TRACE;
        exit(1);
      }
      if(!p_tmp) p_tmp = std::make_shared<Cumulant>(m_edges,m_name+"_tmp",p_varweights);
      p_tmp->Fill(val,weights,0);
      p_tmp->m_fills = ncount;
    }


    void FinishMCB() {
      if(p_tmp) {
        // for(auto w: p_tmp->m_weight) msg_Out()<<w<<" ";
        // msg_Out()<<"\n";
        // for(auto w: p_tmp->m_weights) msg_Out()<<w<<" ";
        // msg_Out()<<"\n";
        p_tmp->m_weight2s = p_tmp->m_weights*p_tmp->m_weights;
        p_tmp->m_barWeight2s = p_tmp->m_barWeights*p_tmp->m_barWeights;
        // for(auto w: p_tmp->m_weight2) msg_Out()<<w<<" ";
        // msg_Out()<<"\n";
        // for(auto w: p_tmp->m_weight2s) msg_Out()<<w<<" ";
        // msg_Out()<<"\n\n";
        *this += *p_tmp;
        p_tmp->Reset();
      }
    }

    void AddZeroPoint(int ncount) {
      DEBUG_FUNC(ncount);
      m_fills += ncount;
    }

    void AddZeroPointMCB(int ncount) {
      DEBUG_FUNC(ncount);
      if(!p_tmp) p_tmp = std::make_shared<Cumulant>(m_edges,m_name,p_varweights);
      p_tmp->m_fills = ncount;
    }

    void Reset() {
      m_weights *= 0;
      m_weight2s *= 0;
      m_barWeights *= 0;
      m_barWeight2s *= 0;
      m_fills = 0;
    }

    void Finalize() {
      FinishMCB();
      if(m_fills > 0) {
        m_weights /= m_fills;
        m_weight2s /= m_fills;
        m_errors = sqrt((m_weight2s-m_weights*m_weights)/(m_fills-1));
        m_barWeights /= m_fills;
        m_barWeight2s /= m_fills;
        m_barErrors = sqrt((m_barWeight2s-m_barWeights*m_barWeights)/(m_fills-1));
      }
      else {
        m_errors = 0.*m_weights;
        m_barErrors = 0.*m_barWeights;
      }
      switch(m_oform) {
      case 1:
        break;
      case 11:
        m_edges = std::log10(m_edges);
        break;
      case 21:
        m_edges = std::log(m_edges);
        break;
      }
    }
    
    void Restore() {
      if(m_fills > 0) {
        m_weights *= m_fills;
        m_weight2s *= m_fills;
        m_barWeights *= m_fills;
        m_barWeight2s *= m_fills;
      }
      switch(m_oform) {
      case 1:
        break;
      case 11:
        m_edges = std::pow(10.,m_edges);
        break;
      case 21:
        m_edges = std::exp(m_edges);
        break;
      }
    }

    void Output(ATOOLS::Gzip_Stream& ofile, const std::string& title) {
      if(!p_varweights) {
        msg_Debugging()<<"No varweights "<<m_name<<".\n";
        InitVarWeights(nullptr);
        DO_STACK_TRACE;
        exit(1);
      }
      DEBUG_FUNC(p_varweights->size());
      for(const auto& var: *p_varweights) {
        *ofile.stream()<<title<<" "<<var.first<<"\n";
        *ofile.stream()<<"# v Sigma{sumW sumW2 err} barSigma{sumW sumW2 err} NumEntries\n";
        const std::valarray<double>& weights = m_weights[var.second];
        const std::valarray<double>& weight2s = m_weight2s[var.second];
        const std::valarray<double>& errors = m_errors[var.second];
        const std::valarray<double>& barWeights = m_barWeights[var.second];
        const std::valarray<double>& barWeight2s = m_barWeight2s[var.second];
        const std::valarray<double>& barErrors = m_barErrors[var.second];

        for(size_t i=0; i<m_nEdges; i++) {
          *ofile.stream()<<m_edges[i]<<" "
                         <<weights[i]<<" "
                         <<weight2s[i]<<" "
                         <<errors[i]<<" "
                         <<barWeights[i]<<" "
                         <<barWeight2s[i]<<" "
                         <<barErrors[i]<<" "
                         <<m_fills<<"\n";
        }
        *ofile.stream()<<"\n";
      }
    }

    const std::string& Name() {return m_name;}

    void InitVarWeights(std::map<std::string,double>* weightsMap=nullptr) {
      DEBUG_FUNC((weightsMap ? weightsMap->size()+1:1));
      std::valarray<size_t> central(m_nEdges);
      std::iota(std::begin(central), std::end(central), 0);
      p_varweights = std::make_shared<VarWeightMap>();
      (*p_varweights)["CENTRAL"] = central;
      if(!weightsMap) {
        m_nVarWeights = 1;
        return;
      }

      m_nVarWeights = weightsMap->size()+1;
      m_nWeights = m_nEdges*m_nVarWeights;

      std::valarray<double> help = std::move(m_weights);
      m_weights = std::valarray<double>(m_nWeights);
      m_weights[central] = std::move(help);

      help = std::move(m_barWeights);
      m_barWeights = std::valarray<double>(m_nWeights);
      m_barWeights[central] = std::move(help);


      help = std::move(m_weight2s);
      m_weight2s = std::valarray<double>(m_nWeights);
      m_weight2s[central] = std::move(help);

      help = std::move(m_barWeight2s);
      m_barWeight2s = std::valarray<double>(m_nWeights);
      m_barWeight2s[central] = std::move(help);


      std::valarray<size_t> idx = central;
      for(const auto& w: *weightsMap) {
        idx += m_nEdges;
        (*p_varweights)[w.first] = idx;
        m_weights[idx] = m_weights[central];
        m_weight2s[idx] = m_weight2s[central];
        m_barWeights[idx] = m_barWeights[central];
        m_barWeight2s[idx] = m_barWeight2s[central];
      }
    }


  protected:
    std::string m_name;
    std::valarray<double> m_edges;
    size_t m_nEdges;
    std::shared_ptr<VarWeightMap> p_varweights = nullptr;
    size_t m_nVarWeights;
    size_t m_nWeights;

    std::valarray<double> m_weights;
    std::valarray<double> m_weight2s;
    std::valarray<double> m_errors;
    std::valarray<double> m_barWeights;
    std::valarray<double> m_barWeight2s;
    std::valarray<double> m_barErrors;

    int m_fills = 0;
    
    // output format, 
    // 1 ... keep bins as is, 
    // 11 .. convert to log10, 
    // 21 ..convert to natural log
    int m_oform = 1;
  
    Cumulant::Ptr p_tmp = nullptr;


  };
}

#endif
