#ifndef RESUM__Bases__Generic_H
#define RESUM__Bases__Generic_H

#include "Tools/CBasis.H"
#include "Tools/CBasis.C"
#include "Tools/CMetric_Base.H"
#include "Bases/Basis_Automate.C"
#include "Bases/Color_Product.C"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Org/Message.H"
#include "Math/asa007.hpp"

using namespace ATOOLS;
using namespace RESUM;



namespace RESUM {

  class CM_Generic : public CMetric_Base  {
    
  private:

    int m_ng,m_nq,m_naq,m_ntot;    
    
  public:
    
    CM_Generic(const CMetric_Key &args);
    void CalcMetric();
    void CalcIMetric();
    void CalcTs();
    void ReadPermutations(ATOOLS::Cluster_Amplitude *ampl);
    std::string m_rpath;
    std::string m_filename;
    
  };
  
}

DECLARE_GETTER(CM_Generic,"CM_Generic",CMetric_Base,CMetric_Key);

CMetric_Base *ATOOLS::Getter<CMetric_Base,CMetric_Key,CM_Generic>::
operator()(const Parameter_Type &args) const 
{
  Flavour_Vector fl;
  for (size_t i(0);i<args.p_ampl->Legs().size();++i)
    if (args.p_ampl->Leg(i)->Flav().Strong())
      fl.push_back(args.p_ampl->Leg(i)->Flav());
  //if (fl.size()<4) return NULL;
  return new CM_Generic(args); 
}

void ATOOLS::Getter<CMetric_Base,CMetric_Key,CM_Generic>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<" CM Generic";
}

CM_Generic::CM_Generic(const CMetric_Key &args):
  CMetric_Base(args), m_ng(0), m_nq(0), m_naq(0), m_ntot(0)
{
  DEBUG_FUNC("");
  for (size_t i(0);i<args.p_ampl->Legs().size();++i) {
    Flavour flav = args.p_ampl->Leg(i)->Flav();
    if (i==0 || i==1) flav=flav.Bar();
    if (flav==Flavour(kf_gluon)) m_ng++;
    if (flav.IsQuark() && !flav.IsAnti()) m_nq++;
    if (flav.IsQuark() && flav.IsAnti()) m_naq++;
  }
  
  m_ntot = m_ng+m_nq+m_naq;
  m_rpath = "../../Bases/pre_calc/";
  m_filename = "";
  for(unsigned i = 0; i < double(m_nq+m_naq)/2.; i++) m_filename = m_filename + "qqb";
  for(unsigned i = 0; i < m_ng; i++) m_filename = m_filename + "g";
  m_map.resize(m_ntot);
  m_pam.resize(m_ntot);
  for (size_t i(0);i<m_ntot;++i) m_pam[m_map[i]=(m_ng+i)%(m_ng+m_naq+m_nq)]=i;
  msg_Debugging()<<"map: "<<m_map<<"\n";
  msg_Debugging()<<"pam: "<<m_pam<<"\n";

  CalcMetric();
  CalcIMetric();
  CalcTs();
  msg_Debugging()<<"Read Permutations...\n";
  ReadPermutations(args.p_ampl);
}

void CM_Generic::ReadPermutations(Cluster_Amplitude *ampl) {
  msg_Debugging()<<*ampl<<"\n";
  int size_connected;
  
  //Get permutations
  msg_Debugging()<<"Read from "<<(m_rpath+m_filename+"/"+m_filename+"_perms.dat").c_str()<<"\n";
  ifstream in( (m_rpath+m_filename+"/"+m_filename+"_perms.dat").c_str() );
  in >> size_connected;
  msg_Debugging()<<"Expect "<<size_connected<<" permutations.\n";
  m_perms.resize(size_connected);
  int temp;
  for (size_t l = 0; l < size_connected; l++){
    msg_Debugging()<<"Perm "<<l<<": ";
    for (size_t m = 0; m < m_ntot; m++){
      in >> temp;
      msg_Debugging()<<ID(ampl->Leg(Map(temp))->Id())<<" ";
      // m_perms[l].push_back(ID(ampl->Leg(Map(temp))->Id()).front());
      m_perms[l].push_back(ID(ampl->Leg(Map(temp))->Id()).front());
    }
    msg_Debugging()<<"\n";
  }
  
  msg_Debugging()<<"Repeat permutations read in:\n";
  for(size_t n = 0; n < m_perms.size(); n++){
    msg_Debugging() << "Perm: " << m_perms[n] << std::endl; 
  }
}

void CM_Generic::CalcMetric(){
 
  int DIM;
  ifstream in( (m_rpath+m_filename+"/"+m_filename+"_met.dat").c_str() );
  in >> DIM;
  
  m_metric.resize(DIM);
  for (size_t i = 0; i < DIM; i++) {
       m_metric[i].resize(DIM);   
       for (size_t j = 0; j < DIM; j++) {
	    in >> m_metric[i][j];
	  }
        }

  in.close();
}


void CM_Generic::CalcIMetric(){
  
  for (size_t i=0;i<m_Imetric.size();i++)
  m_Imetric[i].clear();
  m_Imetric.clear();
  
  m_Imetric = CalcInverse(m_metric);

}

void CM_Generic::CalcTs(){
  
  int DIM;
  ifstream in( (m_rpath+m_filename+"/"+m_filename+".dat").c_str() );
  in >> DIM;
  
  for (size_t i=0;i<m_Tprods.size();i++)
  m_Tprods[i].clear();
  m_Tprods.clear();
  
  for(size_t i = 1; i<=m_ntot; i++){
      for(size_t j = i+1; j<=m_ntot; j++){

     std::vector< std::vector< double > > tmp_tp;
  
     tmp_tp.resize(DIM);
     for (size_t l = 0; l < DIM; l++) {
       tmp_tp[l].resize(DIM);   
       for (size_t m = 0; m < DIM; m++) {
	    in >> tmp_tp[l][m];
	  }
        }
	m_Tprods.push_back(tmp_tp);
      }
    }

}


#endif
