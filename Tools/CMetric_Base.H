#ifndef RESUM__Tools__CMetric_Base_H
#define RESUM__Tools__CMetric_Base_H

#include "Tools/CBasis.H"
#include "PHASIC++/Main/Color_Integrator.H"
#include "ATOOLS/Org/Getter_Function.H"
#include "Math/Matrix.H"
#include <string>

namespace RESUM {
  
  struct CMetric_Key {
  public:
    std::string m_name;
    ATOOLS::Cluster_Amplitude *p_ampl;
    CMetric_Key(const std::string &name,
		ATOOLS::Cluster_Amplitude *const ampl):
      m_name(name), p_ampl(ampl) {}
  };

  class CMetric_Base {
  protected:
    std::string m_name;
    std::vector< std::vector< double > > m_metric,  m_Imetric;
    std::vector<double> m_prefactors;
    std::vector<std::vector<double>> m_trafoMatrix;
    std::vector<  std::vector< std::vector< double > > > m_Tprods;
    std::vector< PHASIC::Idx_Vector > m_perms;
    std::vector< size_t > m_map, m_pam;
  public:
    CMetric_Base();
    CMetric_Base(const CMetric_Key &args);

    void virtual CalcMetric();
    void virtual CalcTs();
    void virtual CalcIMetric();
    
    

    virtual ~CMetric_Base();
    
    static CMetric_Base *GetCM(const CMetric_Key &args);
    
    inline const std::vector<PHASIC::Idx_Vector> &Perms() const { return m_perms; }

    inline size_t Map(const size_t &i) { return m_map[i]; }
    inline size_t Pam(const size_t &i) { return m_pam[i]; }
    
    inline const std::string Name() { return m_name; }

    inline const std::vector<double>& Prefactors()
    { return m_prefactors; }

    MatrixD PrefMatrix() {
      size_t dim = -1;
      if(m_trafoMatrix.empty()) dim = m_metric.size();
      else dim = m_trafoMatrix.at(0).size();
      MatrixD ret(nullptr, dim, dim, 0);
      for(size_t i=0; i<m_prefactors.size(); i++) {
        for(size_t j=0; j<m_prefactors.size(); j++) {
          ret(i,j) = m_prefactors.at(i)*m_prefactors.at(j);
        }
      }
      return ret;
    }

    
    inline const std::vector< std::vector< double > >& CMetric() 
    { return m_metric; }

    inline const std::vector< std::vector< double > >& Imetric() 
    { return m_Imetric; }

    inline const std::vector< std::vector< std::vector< double > > >& Tprods() 
    { return m_Tprods; }

    inline bool hasTrafo() { return !m_trafoMatrix.empty(); }
    
    std::vector<std::vector<double>> TransformationMatrix() {
      if(m_trafoMatrix.empty()) {
        std::vector<std::vector<double>> ret = std::vector<std::vector<double>>(m_metric.size(), std::vector<double>(m_metric.size(),0));
        for(size_t i=0; i<m_trafoMatrix.size(); i++) m_trafoMatrix[i].at(i) = 1;
        return ret;
      }
      return m_trafoMatrix;
    }
  };
  
  typedef std::map<std::string, CMetric_Base*> CMetric_Map;
  typedef std::map<std::string, std::vector< std::vector<double> > > Tprod_Map;
}
#endif
