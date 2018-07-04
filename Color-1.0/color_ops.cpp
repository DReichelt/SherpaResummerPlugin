#ifndef __COLOROPS_CPP
#define __COLOROPS_CPP

#include "color_flow.h"

inline std::vector< Cbasis<int> > CNull(const int& in, const int& ln){
std::vector< Cbasis<int> > results;
Cbasis<int> p1(0,0,1,0);
 
return results;
}


inline std::vector< Cbasis<int> > MC(const int& in, const int& jn, const int& kn,const int& ln){
std::vector< Cbasis<int> > results;
Cbasis<int> p1(2,0,0.5,0);
p1(0,0) = in;
p1(0,1) = ln;
p1(1,0) = kn;
p1(1,1) = jn;

 Cbasis<int> p2(2,0,-0.5,-1);
p2(0,0) = in;
p2(0,1) = jn;
p2(1,0) = kn;
p2(1,1) = ln;

results.push_back(p1); 
results.push_back(p2);

return results;
}


inline std::vector< Cbasis<int> > MPC(const int& io, const int& jo, const int& it, const int& jt,const int& ith, const int& jth){
std::vector< Cbasis<int> > results;
Cbasis<int> p1(3,0,0.5,0);
p1(0,0) = io;
p1(0,1) = jt;
p1(1,0) = ith;
p1(1,1) = jo;
p1(2,0) = it;
p1(2,1) = jth;

Cbasis<int> p2(3,0,-0.5,0);
p2(0,0) = io;
p2(0,1) = jth;
p2(1,0) = it;
p2(1,1) = jo;
p2(2,0) = ith;
p2(2,1) = jt;

results.push_back(p1); 
results.push_back(p2);
 
return results;
}



inline std::vector< Cbasis<int> > PC(const int& io, const int& jo, const int& it, const int& jt,const int& ith, const int& jth,const int& ifo, const int& jfo){
std::vector< Cbasis<int> > results;

Cbasis<int> p1(4,0,0.5,0);
p1(0,0) = io;
p1(0,1) = jt;
p1(1,0) = it;
p1(1,1) = jth;
p1(2,0) = ith;
p1(2,1) = jfo;
p1(3,0) = ifo;
p1(3,1) = jo;

Cbasis<int> p2(4,0,0.5,0);
p2(0,0) = io;
p2(0,1) = jfo;
p2(1,0) = it;
p2(1,1) = jo;
p2(2,0) = ith;
p2(2,1) = jt;
p2(3,0) = ifo;
p2(3,1) = jth;

Cbasis<int> p3(4,0,-0.5,0);
p3(0,0) = io;
p3(0,1) = jth;
p3(1,0) = it;
p3(1,1) = jo;
p3(2,0) = ith;
p3(2,1) = jfo;
p3(3,0) = ifo;
p3(3,1) = jt;

Cbasis<int> p4(4,0,-0.5,0);
p4(0,0) = io;
p4(0,1) = jt;
p4(1,0) = it;
p4(1,1) = jfo;
p4(2,0) = ith;
p4(2,1) = jo;
p4(3,0) = ifo;
p4(3,1) = jth;

results.push_back(p1); 
results.push_back(p2); 
results.push_back(p3); 
results.push_back(p4); 

return results;
}

#endif
