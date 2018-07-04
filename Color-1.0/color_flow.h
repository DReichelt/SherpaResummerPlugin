#ifndef __COLORFLOW_H
#define __COLORFLOW_H

#include <vector>

template <typename T> class Cbasis{
  public:
  unsigned dim;
  unsigned col;
  double constant;
  double pownc;
  std::vector<std::vector<T> > elms;
  Cbasis(unsigned _dim, const T& _initial, double constant, double pownc);
  unsigned _col;
        Cbasis(const Cbasis<T>& rhs);
	virtual ~Cbasis();

  // Operator overloading, for "standard" mathematical matrix operations        
  Cbasis<T>& operator=(const Cbasis<T>& rhs);

  // Transfer operator         
  Cbasis<T> trans(const std::vector<std::vector<T> >& rhs);
  
  // Conjugate object indices
  Cbasis<T> cjgate(const Cbasis<T>& op);

  // Access the individual elements                                             
  T& operator()(const unsigned& dim, const unsigned& col);
  const T& operator()(const unsigned& dim, const unsigned& col) const;

  // glue together two color strucs
  Cbasis<T> cadd(const Cbasis<T>& rhs);

  //Color operator products
  std::vector< std::vector< double > > 
    colPro( std::vector< std::vector< Cbasis<int> > > c, 
  	    std::vector< Cbasis<int> > op);

  void fill4(const int& i, const int& j, const int& k,const int& l);

  void fill6(const int& i, const int& j, const int& k,const int& l, const int& m, const int& n);
  
  void fill8(const int& i, const int& j, const int& k, 
	     const int& l, const int& m, const int& n, const int& o,const int& p);

  void app(const double& newCon, const double& npnc);

  // Get size                                             
  unsigned get_dim() const;
   
  // Get constant                                             
  double get_constant() const;
  
  // Get size                                             
  double get_pownc() const;

  // Print Cbasis object
  void printBasInfo( );
  
};


#include "color_flow.cpp"
#include "color_ops.cpp"
#include "bases.cpp"

#endif
