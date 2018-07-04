#ifndef __BASES_CPP
#define __BASES_CPP

#include "color_flow.h"


inline std::vector< std::vector< Cbasis<int> > > 
Cqqbqqb(const int& in, const int& jn, const int& kn,const int& ln){
std::vector< std::vector< Cbasis<int> > > results;
std::vector< Cbasis<int> > c1,c2; 

Cbasis<int> p(2,0,1,-1);
p.fill4(in,ln,kn,jn); c1.push_back(p);
p.fill4(in,jn,kn,ln); c2.push_back(p);

results.push_back(c1); 
results.push_back(c2);

return results;
}


inline std::vector< std::vector< Cbasis<int> > > 
Cqqbgg(const int& iq, const int& jq, const int& io, const int& jo, const int& it, const int& jt){
std::vector< std::vector< Cbasis<int> > > results;
std::vector< Cbasis<int> > c1,c2,c3; 

double tnorm = sqrt(s_Nc)/(s_Nc*s_Nc-1); 
double dnorm = 1/sqrt(s_Nc*(s_Nc*s_Nc-1));

Cbasis<int> p1(3,0,dnorm,0);
Cbasis<int> p2(3,0,-dnorm,-1);
Cbasis<int> p3(3,0,tnorm,0);
Cbasis<int> p4(3,0,tnorm,0);
Cbasis<int> p5(3,0,-tnorm,-1);
Cbasis<int> p6(3,0,-tnorm,-1);

p1.fill6(iq,jq,io,jt,it,jo);
p2.fill6(iq,jq,io,jo,it,jt); 
p3.fill6(iq,jo,io,jt,it,jq); 
p4.fill6(iq,jt,it,jo,io,jq);
p5.fill6(iq,jo,io,jq,it,jt); 
p6.fill6(iq,jt,it,jq,io,jo); 

c1.push_back(p1); c1.push_back(p2);
c2.push_back(p3); c2.push_back(p5); c2.push_back(p6);
p2.app(tnorm,-2); c2.push_back(p2);

c3.push_back(p4); c3.push_back(p5); c3.push_back(p6); c3.push_back(p2);

results.push_back(c1); 
results.push_back(c2);
results.push_back(c3);

return results;
}


inline std::vector< std::vector< Cbasis<int> > > 
Cgggg(const int& io, const int& jo, const int& it, const int& jt, 
      const int& ith, const int& jth, const int& ifo, const int& jfo){
std::vector< std::vector< Cbasis<int> > > results;
 std::vector< Cbasis<int> > c1, c2, c3, c4, c5, c6, u1pa, u2pa, u3pa; 

 double norm  = 1/sqrt( 2*(8 - 6/(s_Nc*s_Nc) - 3*s_Nc*s_Nc + pow(s_Nc,4) ) );
 double normT = 1/( s_Nc*s_Nc - 1 ); 

Cbasis<int> u1p(4,0,-norm,-1);
u1p.fill8(io,jo,it,jth,ith,jfo,ifo,jt); u1pa.push_back(u1p);
u1p.fill8(io,jo,it,jfo,ifo,jth,ith,jt); u1pa.push_back(u1p);
u1p.fill8(it,jt,io,jth,ith,jfo,ifo,jo); u1pa.push_back(u1p);
u1p.fill8(it,jt,io,jfo,ifo,jth,ith,jo); u1pa.push_back(u1p);
u1p.fill8(ith,jth,io,jt,it,jfo,ifo,jo); u1pa.push_back(u1p);
u1p.fill8(ith,jth,io,jfo,ifo,jt,it,jo); u1pa.push_back(u1p);
u1p.fill8(ifo,jfo,io,jt,it,jth,ith,jo); u1pa.push_back(u1p);
u1p.fill8(ifo,jfo,io,jth,ith,jt,it,jo); u1pa.push_back(u1p);

c1.insert(c1.end(),u1pa.begin(),u1pa.end());
c2.insert(c2.end(),u1pa.begin(),u1pa.end());
c3.insert(c3.end(),u1pa.begin(),u1pa.end());

Cbasis<int> u2p(4,0,2*norm,-2);
u2p.fill8(io,jo,it,jt,ith,jfo,ifo,jth); u2pa.push_back(u2p);
u2p.fill8(io,jo,ith,jth,it,jfo,ifo,jt); u2pa.push_back(u2p);
u2p.fill8(io,jo,ifo,jfo,it,jth,ith,jt); u2pa.push_back(u2p);
u2p.fill8(it,jt,ith,jth,io,jfo,ifo,jo); u2pa.push_back(u2p);
u2p.fill8(it,jt,ifo,jfo,io,jth,ith,jo); u2pa.push_back(u2p);
u2p.fill8(ith,jth,ifo,jfo,io,jt,it,jo); u2pa.push_back(u2p);

c1.insert(c1.end(),u2pa.begin(),u2pa.end());
c2.insert(c2.end(),u2pa.begin(),u2pa.end());
c3.insert(c3.end(),u2pa.begin(),u2pa.end());

Cbasis<int> u3p(4,0,-6*norm,-3);
u3p.fill8(io,jo,it,jt,ith,jth,ifo,jfo);

c1.push_back(u3p);
c2.push_back(u3p);
c3.push_back(u3p);

Cbasis<int> p1(4,0,norm,0);
p1.fill8(io,jt,it,jth,ith,jfo,ifo,jo); c1.push_back(p1);
p1.fill8(io,jfo,ifo,jth,ith,jt,it,jo); c1.push_back(p1);

Cbasis<int> p2(4,0,norm,0);
p2.fill8(io,jt,it,jfo,ifo,jth,ith,jo); c2.push_back(p2);
p2.fill8(io,jth,ith,jfo,ifo,jt,it,jo); c2.push_back(p2);

Cbasis<int> p3(4,0,norm,0);
p3.fill8(io,jth,ith,jt,it,jfo,ifo,jo); c3.push_back(p3);
p3.fill8(io,jfo,ifo,jt,it,jth,ith,jo); c3.push_back(p3);

Cbasis<int> p41(4,0,normT,0);
p41.fill8(io,jth,ith,jo,it,jfo,ifo,jt); c4.push_back(p41);
p41.fill8(io,jfo,ifo,jo,it,jth,ith,jt); c5.push_back(p41);
p41.fill8(io,jt,it,jo,ith,jfo,ifo,jth); c6.push_back(p41);

Cbasis<int> p42(4,0,-normT,-1);
p42.fill8(io,jth,ith,jo,it,jt,ifo,jfo); c4.push_back(p42);
p42.fill8(io,jo,ith,jth,it,jfo,ifo,jt); c4.push_back(p42);
p42.fill8(io,jfo,ifo,jo,it,jt,ith,jth); c5.push_back(p42);
p42.fill8(io,jo,ifo,jfo,it,jth,ith,jt); c5.push_back(p42);
p42.fill8(io,jt,it,jo,ifo,jfo,ith,jth); c6.push_back(p42);
p42.fill8(io,jo,it,jt,ith,jfo,ifo,jth); c6.push_back(p42);

Cbasis<int> p43(4,0,normT,-2);
p43.fill8(io,jo,it,jt,ith,jth,ifo,jfo);
c4.push_back(p43);
c5.push_back(p43);
c6.push_back(p43);

results.push_back(c1); 
results.push_back(c2);
results.push_back(c3);
results.push_back(c4); 
results.push_back(c5);
results.push_back(c6);

return results;
}


#endif
