#ifndef RESUM__Tools__CMetric_Base_H
#define RESUM__Tools__CMetric_Base_H

#include "Tools/CBasis.H"
#include "PHASIC++/Main/Color_Integrator.H"
#include "ATOOLS/Org/Getter_Function.H"
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
    
    inline const std::vector< std::vector< double > > CMetric() 
    { return m_metric; }

    inline const std::vector< std::vector< double > > Imetric() 
    { return m_Imetric; }

    inline const std::vector< std::vector< std::vector< double > > > Tprods() 
    { return m_Tprods; }
  
    
};
  
  typedef std::map<std::string, CMetric_Base*> CMetric_Map;
  typedef std::map<std::string, std::vector< std::vector<double> > > Tprod_Map;
 
}

#endif


