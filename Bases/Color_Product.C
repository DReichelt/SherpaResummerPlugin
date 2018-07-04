#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstdio>
#include <vector>
#include <iostream>
#include "Tools/CBasis.H"
#include "ATOOLS/Phys/Color.H"
#include "ATOOLS/Math/Matrix.H"
using namespace RESUM;
 
// Assign trick
CBasis<int> cjgate(const CBasis<int>& cconj, const CBasis<int>& op) {
  CBasis<int> results(1,0);
  results = cconj;
  
  //Conjugate the fundemental indices
  for( unsigned i = 0; i<cconj.get_Ddim(); i++){
    results.Din[i][0] = fabs(cconj.Din[i][1]);
    results.Din[i][1] = fabs(cconj.Din[i][0]);
  }
   
  //Conjugate the T fund indices
  for( unsigned i = 0; i<cconj.get_Tdim(); i++){
    results.Tin[i][1] = fabs(cconj.Tin[i][2]);
    results.Tin[i][2] = fabs(cconj.Tin[i][1]);
  }
   
  //assign fun-antifundamental indices T <-> F/T
  for (unsigned i=0; i<op.get_Ddim(); i++) {
    for (unsigned j=0; j<cconj.get_Tdim(); j++){
      if( op.Din[i][1]-100 == results.Tin[j][1] ) results.Tin[j][1] = results.Tin[j][1]+100;
      if( op.Din[i][0]-100 == results.Tin[j][2] ) results.Tin[j][2] = results.Tin[j][2]+100;
    }
    for (unsigned j=0; j<cconj.get_Ddim(); j++){
      if( op.Din[i][0]-100 == results.Din[j][1] ) results.Din[j][1] = results.Din[j][1]+100;
      if( op.Din[i][1]-100 == results.Din[j][0] ) results.Din[j][0] = results.Din[j][0]+100;}
  }
 
  
  //assign fun-antifundamental indices T <-> F/T
  for (unsigned i=0; i<op.get_Tdim(); i++) {
    for (unsigned j=0; j<cconj.get_Tdim(); j++){
      if( op.Tin[i][2]-100 == results.Tin[j][1] && results.Tin[j][1]>50) results.Tin[j][1] = results.Tin[j][1]+100;
      if( op.Tin[i][1]-100 == results.Tin[j][2] && results.Tin[j][2]>50) results.Tin[j][2] = results.Tin[j][2]+100;
    }
    for (unsigned j=0; j<cconj.get_Ddim(); j++){
      if( op.Tin[i][1]-100 == results.Din[j][1] ) results.Din[j][1] = results.Din[j][1]+100;
      if( op.Tin[i][2]-100 == results.Din[j][0] ) results.Din[j][0] = results.Din[j][0]+100;}
  }
 
  //assign correct adjoint indices T <-> F/T
  for (unsigned i=0; i<op.get_Tdim(); i++) {
    for (unsigned j=0; j<cconj.get_Tdim(); j++) {  
      if( op.Tin[i][0]-100 == results.Tin[j][0] ) results.Tin[j][0] = results.Tin[j][0]+100;}
    for (unsigned k=0; k<cconj.get_Fdim(); k++) {  
      for(unsigned adj = 0; adj < 3; adj++){
      if( op.Tin[i][0]-100 == results.Fin[k][adj] ) results.Fin[k][adj] = results.Fin[k][adj]+100;}
      }
  }
 
  //avoid double counting internal f-indices
  for (unsigned k=0; k<cconj.get_Fdim(); k++){
    for(unsigned adj = 0; adj < 3; adj++){
      if( results.Fin[k][adj] > 500 ) results.Fin[k][adj] = results.Fin[k][adj]+100;
    }
  }

  //assign indices F <-> F
  for (unsigned i=0; i<op.get_Fdim(); i++){
    for (unsigned j=0; j<cconj.get_Fdim(); j++) {
      for(unsigned adj = 0; adj < 3; adj++){
	if( op.Fin[i][adj]-100 == results.Fin[j][0]) results.Fin[j][0] = results.Fin[j][0]+100;
	if( op.Fin[i][adj]-100 == results.Fin[j][1]) results.Fin[j][1] = results.Fin[j][1]+100;
	if( op.Fin[i][adj]-100 == results.Fin[j][2]) results.Fin[j][2] = results.Fin[j][2]+100;
      }
    }
    for (unsigned j=0; j<cconj.get_Tdim(); j++) {
      for(unsigned adj = 0; adj < 3; adj++){
	if( op.Fin[i][adj]-100 == results.Tin[j][0]) results.Tin[j][0] = results.Tin[j][0]+100;
      }
    }
  }

  /*CBasis<int> tmp(1,0);
  tmp = results;
  Switch order of the T^Adjoint indices to conjugate
  for (unsigned k=0; k<cconj.get_Tdim(); k++){
      results.Tin[k][0] = tmp.Tin[cconj.get_Tdim()-1-k][0];
      }*/
  
  
  return results;
}
 
 
inline std::vector< CBasis<int> >
PC(const int& A, const int& B, const int& D, const int& E)
{
  std::vector< CBasis<int> > results;
  CBasis<int> p1(1,0);
  p1.Fadd(A,B,2001);
  p1.Fadd(2001,D,E);
   
  results.push_back(p1);
  return results;
}
 
 
//identity operator
inline std::vector< CBasis<int> > Ciden(const int& in, const int& ln){
std::vector< CBasis<int> > results;
CBasis<int> p1(1,0);
results.push_back(p1);
return results;
}
 
 
inline std::vector< CBasis<int> > MC(const int& in, const int& jn, const int& kn,const int& ln){
  std::vector< CBasis<int> > results;
  CBasis<int> p1(.5,0);
  p1.Dadd(in,ln);
  p1.Dadd(kn,jn);
 
  CBasis<int> p2(-.5,-1);
  p2.Dadd(in,jn);
  p2.Dadd(kn,ln);
 
  results.push_back(p1);
  results.push_back(p2);
 
  return results;
}
 
inline std::vector< CBasis<int> > MPC(const int& io, const int& jo, const int& it, const int& jt){
std::vector< CBasis<int> > results;
  CBasis<int> p1(1,0);
  p1.Tadd(io,it,25);
  p1.Tadd(jo,25,jt);
 
  CBasis<int> p2(-1,0);
  p2.Tadd(jo,it,25);
  p2.Tadd(io,25,jt);
 
  results.push_back(p1);
  results.push_back(p2);
  
return results;
}
 
inline std::vector< CBasis<int> > PCT(const int& io, const int& it, const int& ith, const int& ifo){
std::vector< CBasis<int> > results;
CBasis<int> p1(2.,0);
 p1.Tadd(io,io,it);
 p1.Tadd(it,it,ith);
 p1.Tadd(ith,ith,ifo);
 p1.Tadd(ifo,ifo,io);
 
CBasis<int> p2(2.,0);
 p2.Tadd(io,io,ifo);
 p2.Tadd(it,it,io);
 p2.Tadd(ith,ith,it);
 p2.Tadd(ifo,ifo,ith);
 
CBasis<int> p3(-2.,0);
 p3.Tadd(io,io,ith);
 p3.Tadd(it,it,io);
 p3.Tadd(ith,ith,ifo);
 p3.Tadd(ifo,ifo,it);
 
CBasis<int> p4(-2.,0);
 p4.Tadd(io,io,it);
 p4.Tadd(it,it,ifo);
 p4.Tadd(ith,ith,io);
 p4.Tadd(ifo,ifo,ith);
 
results.push_back(p1);
results.push_back(p2);
results.push_back(p3);
results.push_back(p4);
 
return results;
}
 
 
//write to expression T's
void Texpre(CBasis<int> ct, ATOOLS::Expression &ex, int &pl){
  int pl2 = 0;
  int tA, tf, taf;
  ex.push_back(ATOOLS::CNumber::New(ct.get_constant()*pow(s_Nc,ct.get_pownc())));
   
  for(unsigned i = 0; i < ct.get_Tdim(); i++){
    tA = ct.Tin[i][0];
    tf  = ((fabs(ct.Tin[i][1]-50) > 10 && fabs(ct.Tin[i][1]-150) > 10) ? 1000+pl : 0) + fabs(ct.Tin[i][1]);
    taf = ((fabs(ct.Tin[i][2]-50) > 10 && fabs(ct.Tin[i][2]-150) > 10) ? 1000+pl : 0) + fabs(ct.Tin[i][2]);
    ex.push_back(ATOOLS::Fundamental::New(tA,tf,taf));  
    pl2++;
  }
 
  for(unsigned i = 0; i < ct.get_Ddim(); i++){
    ex.push_back(ATOOLS::Delta::New(fabs(ct.Din[i][0]),fabs(ct.Din[i][1])));
  }
 
  for(unsigned i = 0; i < ct.get_Fdim(); i++){
    ex.push_back(ATOOLS::Adjoint::New(ct.Fin[i][0],ct.Fin[i][1],ct.Fin[i][2]));
  }
 
  pl = 2*(pl+pl2)+10;
}
 

// Normalize the basis
void normalize_basis(ABasis &unc){
 
  double hold,normfac;
   
  for(unsigned i = 0; i<unc.size(); i++){
      normfac=0;
      for(unsigned k = 0; k<unc[i].size(); k++){
        for(unsigned l = 0; l<unc[i].size(); l++){
        ATOOLS::Expression expression;
        expression.clear();
        int place = 0;
         
        Texpre(unc[i][k], expression, place);
        Texpre(cjgate(unc[i][l],Ciden(0,0)[0]),expression, place);
 
	expression.Evaluate();
 
        //append the i-th element
        hold = real(expression.Result());
        normfac = normfac + hold;
        }
            }
            for(unsigned m = 0; m < unc[i].size(); m++){
        unc[i][m].app(unc[i][m].get_constant()*(1/sqrt(normfac)),unc[i][m].get_pownc());
        }
       }
  return;
}
 
 
// evaluates color product <c|T_i T_j|c> for either
// (set to fundamental here)
 
std::vector< std::vector< double > > colProT( ABasis c , std::vector< CBasis<int> > op, double v_NC = 3.0 )
{
  double hold = 0;
  unsigned test = 0;
  CBasis<int> take(1,0);
  std::vector< std::vector<double> > mat( c.size() );
   
  for(unsigned i = 0; i<c.size(); i++){
    mat[i].resize(c.size());
    for(unsigned j = i; j<c.size(); j++){
      mat[j].resize(c.size());
      for(unsigned k = 0; k<c[i].size(); k++){
        for(unsigned l = 0; l<c[j].size(); l++){
        for(unsigned m = 0; m<op.size(); m++){
        ATOOLS::Expression expression;
        expression.clear();
	expression.SetNC(v_NC);
        int place = 0;
         
        Texpre(c[i][k], expression, place);
        Texpre(op[m], expression, place);
        Texpre(cjgate(c[j][l],op[m]), expression, place);
       
       	expression.Evaluate();
	
        //build Nc = 3 matrix
        hold = real(expression.Result()) + imag(expression.Result());
        mat[i][j] += hold;
        if( fabs( mat[i][j] ) < s_eps ){hold=0; mat[i][j] = 0; mat[j][i] = 0; }
             if(i!=j) mat[j][i] += hold;
       }
       }
      }
     }
    }
  return mat;
}
