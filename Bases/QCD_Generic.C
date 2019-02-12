#ifndef RESUM__Bases__Generic_H
#define RESUM__Bases__Generic_H

#include "Tools/CBasis.H"
#include "Tools/CBasis.C"
#include "Tools/CMetric_Base.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Org/Message.H"
#include "Math/asa007.hpp"
#include "ATOOLS/Org/Run_Parameter.H"

//#include "ATOOLS/Org/Data_Reader.H"


#include "Tools/Reader.H"
#include "Tools/Files.H"

#include <algorithm>
#include <regex>
#include <locale>

using namespace ATOOLS;
using namespace RESUM;



namespace RESUM {

  template <typename T>
    inline std::vector<std::vector<T>> VectorToMatrix(const std::vector<T>& vec, size_t r, size_t c) {
    if(vec.size() != r*c) THROW(fatal_error, "Wrong dimension for vector.");
    std::vector<std::vector<double>> ret(r);
    for(size_t row=0; row<ret.size(); row++) {
      ret[row] = {vec.begin()+c*row, vec.begin()+c*(row+1)};
    }
    return ret;
  }

  
  template <typename T>
  inline std::vector<std::vector<T>> VectorToMatrix(const std::vector<T>& vec, size_t n) {
    /* if(vec.size() != n*n) THROW(fatal_error, "Wrong dimension for vector."); */
    /* std::vector<std::vector<double>> ret(n); */
    /* for(size_t row=0; row<ret.size(); row++) { */
    /*   ret[row] = {vec.begin()+n*row, vec.begin()+n*(row+1)}; */
    /* } */
    return VectorToMatrix(vec,n,n);
  }


  class CM_Generic : public CMetric_Base  {
    
  private:
    Reader m_reader;
    int m_ng,m_nq,m_naq,m_ntot;    
    
  public:
    
    CM_Generic(const CMetric_Key &args);
    void CalcMetric();
    void CalcIMetric();
    void CalcTs();
    void CalcTransformationMatrix();
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
  m_filename = "";
  for(unsigned i = 0; 2*i < m_nq+m_naq; i++) {
    m_filename = m_filename + "qqb";
  }
  for(unsigned i = 0; i < m_ng; i++) m_filename = m_filename + "g";
  m_map.resize(m_ntot);
  m_pam.resize(m_ntot);
  for (size_t i(0);i<m_ntot;++i) m_pam[m_map[i]=(m_ng+i)%(m_ng+m_naq+m_nq)]=i;

  // TODO: look in more reasonable localtions and terminate properly if
  //       nothing is found
  m_rpath = "./"+rpa->gen.Variable("RESUM::pre_calc")+"/";
  ifstream test1((m_rpath+m_filename+".dat").c_str());
  if(!test1.good()){
    msg_Debugging()<<"No file "<<m_filename<<".dat found in "<<m_rpath<<".\n";
    m_rpath = RESUM::FILENAMES::SHARE_DIR + "/" + rpa->gen.Variable("RESUM::pre_calc")+"/";
  }
  ifstream test2((m_rpath+m_filename+".dat").c_str());
  if(!test2.good()){
    THROW(fatal_error, "Did not find file "+m_filename+".dat in "+m_rpath);
  }

  msg_Debugging()<<"File "<<m_filename<<".dat found in "<<m_rpath<<".\n";

  m_reader.comments = {"#","%","//"};
  m_reader.file_path = {m_rpath+"/"+m_filename+".dat"};
  m_reader.read();

  msg_Debugging()<<"map: "<<m_map<<"\n";
  msg_Debugging()<<"pam: "<<m_pam<<"\n";

  CalcMetric();
  CalcIMetric();
  CalcTs();
  CalcTransformationMatrix();
  msg_Debugging()<<"Read Permutations...\n";
  ReadPermutations(args.p_ampl);
}


void CM_Generic::ReadPermutations(Cluster_Amplitude *ampl) {
  msg_Debugging()<<*ampl<<"\n";

  
  //Get permutations
  msg_Debugging()<<"Read from "<<(m_rpath+m_filename+"/"+m_filename+"_perms.dat").c_str()<<"\n";


  int size_connected = m_reader.GetValue<int>("SIZE_CONNECTED",-1);
  msg_Debugging()<<size_connected<<"\n";

  if (size_connected < 0) {
    msg_Debugging()<<"Still old format.\n";
    ifstream in( (m_rpath+m_filename+"/"+m_filename+"_perms.dat").c_str() );
    in >> size_connected;
    
    msg_Debugging()<<"Expect "<<size_connected<<" permutations.\n";
    m_perms.clear();
    m_perms.reserve(size_connected);
    std::string line;
    size_t l = 0;
    while(std::getline(in,line)) {
      if(l >= size_connected) break;
      if(line.empty()) continue;
      l++;
      size_t begin = line.find("(");
      size_t end = line.find(")");
      msg_Debugging()<<"Read Line: "<<line<<", perm starts at "<<begin<<" and ends at "<<end<<"\n";
      if(begin != 0 && begin != std::string::npos) {

      std::string prefactor = line.substr(0,begin-1);
      //prefactor.erase(std::remove_if(prefactor.begin(), prefactor.end(), std::isspace), prefactor.end());
      msg_Debugging()<<"Add prefactor "<<prefactor<<"\n";
      m_prefactors.push_back(std::stod(prefactor));
      }
      else {
        m_prefactors.push_back(1);
      }
      if(begin == std::string::npos) begin=0;
      else begin += 1;
      if(end == std::string::npos) end=line.size()+1;
      msg_Debugging()<<"Perm "<<line.substr(begin,end-begin)<<".\n";
      std::vector<std::string> tmp = split(line.substr(begin,end-begin)," ");
      m_perms.emplace_back(PHASIC::Idx_Vector());
      m_perms.back().reserve(m_ntot);
      msg_Debugging()<<"Add permutation "<<l<<": ";
      for(std::string t: tmp) {
        msg_Debugging()<<std::stoi(t)<<" -> "<<Map(std::stoi(t))<<" -> "<<ID(ampl->Leg(Map(std::stoi(t)))->Id()).front()<<"; ";
        m_perms.back().emplace_back(ID(ampl->Leg(Map(std::stoi(t)))->Id()).front());
      }
      msg_Debugging()<<"\n";
    }
  }
  else {
    for(int i=0; i<size_connected; i++) {
      string pref = "a_"+std::to_string(i);
      string amp = "A_"+std::to_string(i);
      double prefactor = m_reader.GetValue<double>(pref,1);
      msg_Debugging()<<"Add prefactor "<<pref<<" = "<<prefactor<<"\n";
      m_prefactors.push_back(prefactor);
      std::vector<int> tmp = m_reader.GetVector<int>(amp);
      msg_Debugging()<<"Add permutation "<<amp<<" = "<<tmp<<"\n";
      m_perms.emplace_back(PHASIC::Idx_Vector());
      for(const int t: tmp) {
        msg_Debugging()<<t<<" -> "<<Map(t)<<" -> "<<ID(ampl->Leg(Map(t))->Id()).front()<<"; ";
        m_perms.back().emplace_back(ID(ampl->Leg(Map(t))->Id()).front());
      }
      msg_Debugging()<<"\n";
    }
  }
    
  
  // int temp;
  // for (size_t l = 0; l < size_connected; l++){
  //   msg_Debugging()<<"Perm "<<l<<": ";
  //   for (size_t m = 0; m < m_ntot; m++){
  //     in >> temp;
  //     msg_Debugging()<<temp<<" -> "<<Map(temp)<<" -> "<<ID(ampl->Leg(Map(temp))->Id())<<" ";
  //     m_perms[l].push_back(ID(ampl->Leg(Map(temp))->Id()).front());
  //   }
  //   msg_Debugging()<<"\n";
  // }
  
  msg_Debugging()<<"Repeat permutations read in:\n";
  for(size_t n = 0; n < m_perms.size(); n++){
    msg_Debugging() <<"Prefactor = "<<m_prefactors.at(n)<< ", Perm = " << m_perms.at(n) << std::endl; 
  }
}

void CM_Generic::CalcMetric(){
  int DIM = m_reader.GetValue<int>("DIM_CS",-1);
  msg_Debugging()<<DIM<<"\n";
  if(DIM<0) DIM = m_reader.GetValue<int>("DIM",-1);

  if(DIM < 0) {
    THROW(fatal_error, "Old format not supported anymore.")
    /* msg_Debugging()<<"Read old format from "<<m_rpath+m_filename+'/'+m_filename+"_met.dat"<<"\n"; */
    /* ifstream in( (m_rpath+m_filename+'/'+m_filename+"_met.dat").c_str() ); */
    /* in >> DIM; */
    
    /* m_metric.resize(DIM); */
    /* for (size_t i = 0; i < DIM; i++) { */
    /*   m_metric[i].resize(DIM);    */
    /*   for (size_t j = 0; j < DIM; j++) { */
    /*     in >> m_metric[i][j]; */
    /*   } */
    /* } */
    /* in.close(); */
  }
  else {
    msg_Debugging()<<"Read in metric.\n";
    std::vector<double> met = m_reader.GetVector<double>("METRIC");
    m_metric = VectorToMatrix(met, DIM);
  }
}


void CM_Generic::CalcIMetric(){

  m_Imetric = InverseLowerTriangular(m_metric);
  /* for (size_t i=0;i<m_Imetric.size();i++) */
  /* m_Imetric[i].clear(); */
  /* m_Imetric.clear(); */
  
  /* m_Imetric = CalcInverse(m_metric); */

}

void CM_Generic::CalcTs(){
  
  int DIM = m_reader.GetValue<int>("DIM_CS", -1);
  if(DIM < 0) DIM = m_reader.GetValue<int>("DIM", -1);
  if(DIM < 0) {
    ifstream in( (m_rpath+m_filename+"/"+m_filename+".dat").c_str() );
    in >> DIM;
    
    for (size_t i=0;i<m_Tprods.size();i++)
    m_Tprods[i].clear();
    m_Tprods.clear();
    
    for(int i = 1; i<=m_ntot; i++){
      for(int j = i+1; j<=m_ntot; j++){
        
        std::vector< std::vector< double > > tmp_tp;
        
        tmp_tp.resize(DIM);
        for (int l = 0; l < DIM; l++) {
          tmp_tp[l].resize(DIM);   
          for (int m = 0; m < DIM; m++) {
	    in >> tmp_tp[l][m];
	  }
        }
	m_Tprods.push_back(tmp_tp);
      }
    }
  }
  else {
    m_Tprods.clear();
    for(int i=0; i<m_ntot; i++) {
      for(int j=i+1; j<m_ntot; j++) {
        string mat = "C_"+std::to_string(i)+std::to_string(j);
        msg_Debugging()<<"Read matrix "<<mat<<".\n";
        m_Tprods.push_back(m_reader.GetMatrix<double>(mat,DIM));
        if(m_Tprods.back().empty() || m_Tprods.back().size()!=DIM ||
           m_Tprods.back().back().size()!=DIM) {
          // if not all Tprods are found, or they are not complete,
          // it doesnt make sense to have some of
          // them, and code that uses them can not proceed
          m_Tprods.clear();
          return;
        }
      }
    }
  }
}

void CM_Generic::CalcTransformationMatrix() {
  int DIM_BASIS = m_reader.GetValue<int>("DIM_BASIS", -1);
  int DIM_CS = m_reader.GetValue<int>("DIM_CS", -1);
  m_trafoMatrix.clear();
  if(!(DIM_BASIS < 0 || DIM_CS < 0)) {
    m_trafoMatrix = m_reader.GetMatrix<double>("TRAFO",DIM_BASIS); 
    if(m_trafoMatrix.size()!=DIM_CS || m_trafoMatrix.back().size()!=DIM_BASIS) {
      // if the matrix was not complete, treat as if none was found
      m_trafoMatrix.clear();
    }
  }
}

#endif
