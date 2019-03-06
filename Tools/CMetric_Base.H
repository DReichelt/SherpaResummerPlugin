#ifndef RESUM__Tools__CMetric_Base_H
#define RESUM__Tools__CMetric_Base_H

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
    
    MatrixD m_metric,  m_Imetric;

    std::vector<double> m_prefactors;
    
    MatrixD m_trafoMatrix;
    std::vector<MatrixD> m_Tprods;
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
      if(!hasTrafo()) dim = m_metric.numCols();
      else dim = m_trafoMatrix.numRows();
      MatrixD ret(nullptr, dim, dim, 0);
      for(size_t i=0; i<m_prefactors.size(); i++) {
        for(size_t j=0; j<m_prefactors.size(); j++) {
          ret(i,j) = m_prefactors.at(i)*m_prefactors.at(j);
        }
      }
      return ret;
    }
    
    inline const MatrixD& CMetric() const { return m_metric; }

    inline const MatrixD& Imetric() const { return m_Imetric; }
    
    inline const std::vector<MatrixD>& Tprods() const { return m_Tprods; }

    inline bool hasTrafo() const { return m_trafoMatrix.numElements()>0; }
    
    MatrixD TransformationMatrix() {
      if(!hasTrafo()) {
        return MatrixD::diagonal(1,m_metric.numCols(),0);
      }
      return m_trafoMatrix;
    }
  };
  
  typedef std::map<std::string, CMetric_Base*> CMetric_Map;
  //  typedef std::map<std::string, std::vector< std::vector<double> > > Tprod_Map;
}
#endif
