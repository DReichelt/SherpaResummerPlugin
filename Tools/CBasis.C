#ifndef RESUM_CBASIS_C
#define RESUM_CBASIS_C

#include "Tools/CBasis.H"
#include "ATOOLS/Org/Message.H"
#include <iostream>
#include <math.h>

using namespace RESUM;
using namespace std;

const double s_Nc = 3.;
const double s_CF = (s_Nc*s_Nc-1.)/(2.*s_Nc);
const double s_CA = s_Nc;
const double s_TR = 1./2.;
const double s_eps = .000001;


// Parameter Constructor                                                                      
template<typename T>
CBasis<T>::CBasis(double _constant, double _pownc) { 
  constant = _constant;
  pownc = _pownc;
  }

// Copy Constructor                                                                      
template<typename T>
CBasis<T>::CBasis(const CBasis<T>& rhs) {
  constant  = rhs.get_constant();
  pownc  = rhs.get_pownc();
  Tin = rhs.Tin;
  Din = rhs.Din;
  Fin = rhs.Fin;
  }

// (Virtual) Destructor                                                                        
  template<typename T>
  CBasis<T>::~CBasis() {}


// Assignment Operator  
template<typename T>
CBasis<T>& CBasis<T>::operator=(const CBasis<T>& rhs) {
  if (&rhs == this)
  return *this;

  double new_constant = rhs.get_constant();
  double new_pownc    = rhs.get_pownc();
  constant = new_constant;
  pownc = new_pownc;

  unsigned new_T_dim = rhs.get_Tdim();
  unsigned new_D_dim = rhs.get_Ddim();
  unsigned new_F_dim = rhs.get_Fdim();
 
  Tin.resize(new_T_dim);
  for (unsigned i=0; i<Tin.size(); i++) {
    Tin[i].resize(3);
    for (unsigned j=0; j<Tin[i].size(); j++) {
      Tin[i][j] = rhs.Tin[i][j];
    }
  }
    
  Din.resize(new_D_dim);
  for (unsigned i=0; i<Din.size(); i++) {
    Din[i].resize(2);
    for (unsigned j=0; j<Din[i].size(); j++) {
      Din[i][j] = rhs.Din[i][j];
    }
  }
  
  Fin.resize(new_F_dim);
  for (unsigned i=0; i<Fin.size(); i++) {
    Fin[i].resize(3);
    for (unsigned j=0; j<Fin[i].size(); j++) {
      Fin[i][j] = rhs.Fin[i][j];
    }
  }

  return *this;
  }


template<typename T>
  unsigned CBasis<T>::get_Tdim() const {
  return this->Tin.size();}
template<typename T>
  unsigned CBasis<T>::get_Ddim() const {
  return this->Din.size();}
template<typename T>
  unsigned CBasis<T>::get_Fdim() const {
  return this->Fin.size();}

template<typename T>
  double CBasis<T>::get_constant() const {
  return this->constant;
}

template<typename T>
  double CBasis<T>::get_pownc() const {
  return this->pownc;
}


template<typename T>
  CBasis<T>& CBasis<T>::Tadd(const T& adj,const T& fun,const T& afun){
  std::vector<T> tares(3);
  tares[0] = adj;
  tares[1] = fun;
  tares[2] = afun;
  Tin.push_back(tares);
  return *this;
  }

template<typename T>
  CBasis<T>& CBasis<T>::Dadd(const T& fun,const T& afun){
  std::vector<T> tares(2);
  tares[0] = fun;
  tares[1] = afun;
  Din.push_back(tares);
  return *this;
  }

template<typename T>
  CBasis<T>& CBasis<T>::Fadd(const T& adj1,const T& adj2,const T& adj3){
  std::vector<T> tares(3);
  tares[0] = adj1;
  tares[1] = adj2;
  tares[2] = adj3;
  Fin.push_back(tares);
  return *this;
  }



// glue together two color strucs                                   
template<typename T>
CBasis<T> CBasis<T>::cadd( const CBasis<T>& rhs ) {
  CBasis result(constant*rhs.get_constant(), pownc + rhs.get_pownc());
  result.Tin.insert(result.Tin.end(),Tin.begin(),Tin.end());
  result.Tin.insert(result.Tin.end(),rhs.Tin.begin(),rhs.Tin.end());
  result.Din.insert(result.Din.end(),Din.begin(),Din.end());
  result.Din.insert(result.Din.end(),rhs.Din.begin(),rhs.Din.end());
  result.Fin.insert(result.Fin.end(),Fin.begin(),Fin.end());
  result.Fin.insert(result.Fin.end(),rhs.Fin.begin(),rhs.Fin.end());
  return result;
}


//Print a Cbasis element
template <typename T>
void CBasis<T>::printBasInfo() 
{
  msg_Debugging() << constant << "*Nc^" << pownc << " " << std::endl;
    if(Tin.size() > 0){
      msg_Debugging() << "T.fund" << std::endl;
    for (int i=0; i<3; i++) { 
    msg_Debugging() << "[ ";
    for(int j=0; j<Tin.size(); j++){
      if( Tin[j][i] > 0) msg_Debugging() << " ";
      msg_Debugging() << Tin[j][i] << " ";
      }
    msg_Debugging() << "]";
    msg_Debugging() << std::endl;
      }
    }

    if(Din.size() > 0){
    msg_Debugging() << "D.fund" << std::endl;
    for (int i=0; i<2; i++) { 
    msg_Debugging() << "[ ";
    for(int j=0; j<Din.size(); j++){
      if( Din[j][i] > 0) msg_Debugging() << " ";
      msg_Debugging() << Din[j][i] << " ";
      }
    msg_Debugging() << "]";
    msg_Debugging() << std::endl;
      }
    }
    
    if(Fin.size() > 0){
    msg_Debugging() << "F.Adjoint" << std::endl;
    for (int i=0; i<3; i++) { 
    msg_Debugging() << "[ ";
    for(int j=0; j<Fin.size(); j++){
      if( Fin[j][i] > 0) msg_Debugging() << " ";
      msg_Debugging() << Fin[j][i] << " ";
      }
    msg_Debugging() << "]";
    msg_Debugging() << std::endl;
      }
    }
    msg_Debugging() << std::endl;
}


template<typename T>
void CBasis<T>::app(const double& newCon, const double& npnc){
  this->constant=newCon;
  this->pownc=npnc;
}
 

//------------------------------------------------------------
//------------------------------------------------------------
//------------------------------------------------------------

//Non-members

//Print a matrix
template <typename T>
void printMat(std::vector< std::vector <T> > &rhs) {
  msg_Debugging() << std::endl;
  for (unsigned i=0; i<rhs.size(); i++) {
    msg_Debugging() << "{";
    for (unsigned j=0; j<rhs[i].size(); j++) {
      msg_Debugging() << (fabs(rhs[i][j]) > .0000001 ? rhs[i][j] : 0);
      if(j<rhs[j].size()-1) msg_Debugging() << ",";
    }
    msg_Debugging() << "},";
    msg_Debugging() << std::endl;
  }
  msg_Debugging() << std::endl;
}
 
 
// returns minus a matrix
inline void Minus(std::vector< std::vector <double> > &ts){
  for(unsigned i = 0; i<ts.size(); i++){
    for(unsigned j = 0; j<ts[i].size(); j++){
    ts[i][j] = -ts[i][j];
    }
  }
}

#endif


