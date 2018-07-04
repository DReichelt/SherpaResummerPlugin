#ifndef Tools_H
#define Tools_H

#if defined __GNUC__
#if __GNUC__ > 2
#include <sstream>
#include <iomanip>
typedef std::stringstream MyStrStream;
#else
#include <strstream>
#include <string>
#include <iomanip>
typedef std::strstream MyStrStream;
#endif
#endif
#if defined __sgi && !defined __GNUC__
#include <sstream>
#include <iomanip>
typedef std::stringstream MyStrStream;
#endif
#include <sys/times.h>
#include <time.h>

#define THROW(TAG,MESSAGE) \
  { std::cerr<<"error: "<<MESSAGE<<std::endl; exit(1); }

#define msg_Error() std::cerr
#define msg_Info() std::cout
#define msg_Out() std::cout
#define mm_up(NUMBER) "\e["<<NUMBER<<"A"
#define bm_cr "\r"

namespace ATOOLS
{

template <class Value_Type>
std::string ToString(const Value_Type value,
                     const size_t precision=12)
{
    MyStrStream converter;
    std::string converted;
    converter.precision(precision);
    converter<<value;
    converter>>converted;
    return converted;
}

template <class Value_Type>
Value_Type ToType(const std::string &value,
                  const size_t precision=12)
{
    MyStrStream converter;
    Value_Type converted;
    converter.precision(precision);
    converter<<value;
    converter>>converted;
    return converted;
}

inline double UserTime()
{
    tms ctime;
    times(&ctime);
    double clk_tck(100.0);
    return ctime.tms_utime/clk_tck;
}


template <class Type>
const Type &MaxT(const Type &a,const Type &b)
{
    return a>b?a:b;
}

}// end of namespace ATOOLS;

#endif
