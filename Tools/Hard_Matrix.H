#ifndef RESUM__Tools__Hard_Matrix_H
#define RESUM__Tools__Hard_Matrix_H

#include "PHASIC++/Main/Color_Integrator.H"
#include "ATOOLS/Math/MyComplex.H"

#include <iostream>

namespace RESUM {

  struct Hard_Matrix: public std::vector<std::vector<Complex> > 
  {
  public:
    std::vector<PHASIC::Idx_Vector> m_id;
    Hard_Matrix(const size_t &n=0):
    std::vector<std::vector<Complex> >
    (n,std::vector<Complex>(n,Complex(0.0))), m_id(n) {}
  };

  
  std::ostream &operator<<(std::ostream &s,const Hard_Matrix &h);
  int size( const Hard_Matrix &h);
  int bas( const Hard_Matrix &h , const int i,const int k);
  double access(const Hard_Matrix &h,const int i,const int j);
}

#endif

