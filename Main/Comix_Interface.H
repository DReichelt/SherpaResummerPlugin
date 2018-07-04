#ifndef RESUM__Tools__Comix_Interface_H
#define RESUM__Tools__Comix_Interface_H

#include "PHASIC++/Process/Process_Base.H"
#include "METOOLS/Main/Spin_Structure.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "Tools/Hard_Matrix.H"

namespace RESUM {

  struct Hard_Coefficient {
  public:

    PHASIC::Idx_Vector m_id;

    METOOLS::Spin_Amplitudes m_amps;

    inline Hard_Coefficient
    (const PHASIC::Idx_Vector &id,
     const METOOLS::Spin_Amplitudes &amps):
      m_id(id), m_amps(amps) {}

  };

  std::ostream &operator<<(std::ostream &s,const Hard_Coefficient &c);

  class Comix_Interface {
  private:

    Hard_Matrix *p_h;

    PHASIC::NLOTypeStringProcessMap_Map m_pmap;
    PHASIC::Process_Vector m_procs;

    PHASIC::Process_Base *GetProcess
    (ATOOLS::Cluster_Amplitude *const ampl);

  public:

    Comix_Interface();

    ~Comix_Interface();

    Hard_Matrix *ComputeHardMatrix
    (ATOOLS::Cluster_Amplitude *const ampl,
     const std::vector<PHASIC::Idx_Vector> &perms);

    void Reset();

  };// end of class Comix_Interface

}// end of namespace RESUM

#endif
