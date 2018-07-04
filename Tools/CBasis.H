#ifndef RESUM__Tools__CBasis_H
#define RESUM__Tools__CBasis_H

#include <vector>

namespace RESUM {

template <typename T> class CBasis{
  public:
  double constant;
  double pownc;
  std::vector< std::vector<T> > Tin;
  std::vector< std::vector<T> > Din;
  std::vector< std::vector<T> > Fin;
 
  //Constructor
  CBasis(double constant, double pownc);
        CBasis(const CBasis<T>& rhs);
	virtual ~CBasis();
       
  //Assignment operator
  CBasis<T>& operator=(const CBasis<T>& rhs);

  //Access indices
  T& Tacc(const unsigned& dims, const unsigned& cols);
  T& Dacc(const unsigned& dims, const unsigned& cols);
  T& Facc(const unsigned& dims, const unsigned& cols);
  
  //Glue together two basis elements
  CBasis<T> cadd( const CBasis<T>& rhs );
  
  //Add color tensors
  CBasis<T>& Tadd(const T& adj,const T& fun,const T& afun);
  CBasis<T>& Dadd(const T& fun,const T& afun);
  CBasis<T>& Fadd(const T& adj1,const T& adj2,const T& adj3);
  
  //Conjugate basis element
  CBasis<T> cjgate(const CBasis<T>& op);
  
  //Append the constants and power of Nc
  void app(const double& newCon, const double& npnc);
  
  //Append the constants and power of Nc
  unsigned get_Tdim() const;
  unsigned get_Ddim() const;
  unsigned get_Fdim() const;
                                
  double get_constant() const;                                  
  double get_pownc() const;
  
  // Print Cbasis element
  void printBasInfo( );

};


typedef std::vector< std::vector< CBasis<int> > > ABasis;

}


#endif
