#include "Color.H"
#include "Tools.H"
#include <termios.h>
#include <unistd.h>
#ifdef MALLOC_TRACE
#error MALLOC_TRACE should not be defined
#include <mcheck.h>
#endif

using namespace ATOOLS;

int main(int argc,char **argv)
{

#ifdef MALLOC_TRACE
  setenv("MALLOC_TRACE","malloc_trace_color.log",1);
  mtrace();
#endif
  
  std::string expr;
  expr = "f_[1,4,501]*f_[501,3,502]*f_[502,2,5]*t_[3,11,12]*t_[4,13,11]*t_[5,12,13]*t_[1,14,15]*t_[2,15,14]";
  
  for (int i=1;i<argc;++i) expr+=argv[i];
  Expression expression(expr);
  expression.Print();
  expression.Evaluate();
  std::cout<<"Color: calculating -> "<<expr<<" = "
	   <<expression.Result()<<std::endl;
  
  std::cout << std::endl;
  std::cout << "=======================================" << std::endl;
  std::cout << "=======================================" << std::endl;
  std::cout << "=======================================" << std::endl;
  std::cout << std::endl;

  
  ATOOLS::Expression ex;
  ex.clear();
  ex.push_back(ATOOLS::Adjoint::New(1,4,501));
  ex.push_back(ATOOLS::Adjoint::New(501,3,502));
  ex.push_back(ATOOLS::Adjoint::New(502,2,5));

  ex.push_back(ATOOLS::Fundamental::New(3,11,12));
  ex.push_back(ATOOLS::Fundamental::New(5,13,11));
  ex.push_back(ATOOLS::Fundamental::New(4,12,13));
  ex.push_back(ATOOLS::Fundamental::New(1,14,15));
  ex.push_back(ATOOLS::Fundamental::New(2,15,14));

  ex.Print();
  ex.Evaluate();
  std::cout<<"Color: calculating -> "<<" = "
  <<ex.Result()<<std::endl;


  return 0;

#ifdef MALLOC_TRACE
  muntrace();
#endif

}
