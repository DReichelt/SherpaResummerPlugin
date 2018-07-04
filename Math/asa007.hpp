# include <vector>
#include "ATOOLS/Math/Matrix.H"
std::vector< std::vector<double> > CalcInverse( std::vector< std::vector<double> > &cmet);
void cholesky ( double a[], int n, int nn, double u[], int *nullty, 
  int *ifault );
double r8_abs ( double x );
void syminv ( double a[], int n, double c[], double w[], int *nullty, 
  int *ifault );
void timestamp ( void );
