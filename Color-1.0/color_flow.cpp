#ifndef __COLORFLOW_CPP
#define __COLORFLOW_CPP

#include "color_flow.h"

using namespace std;

const double s_Nc = 3.;
const double s_CF = (s_Nc*s_Nc-1.)/(2.*s_Nc);
const double s_CA = s_Nc;
const double s_TR = 1./2.;

// Parameter Constructor                                                                      

template<typename T>
Cbasis<T>::Cbasis(unsigned _dim, const T& _initial, double _constant, double _pownc) {
  elms.resize(_dim);
    for (unsigned i=0; i<elms.size(); i++) {
    elms[i].resize(2, _initial);
    }
  dim = _dim;
  constant = _constant;
  pownc = _pownc;
  }


// Copy Constructor                                                                                                                                      
template<typename T>
Cbasis<T>::Cbasis(const Cbasis<T>& rhs) {
  elms = rhs.elms;
  dim  = rhs.get_dim();
  constant  = rhs.get_constant();
  pownc  = rhs.get_pownc();
  }


// (Virtual) Destructor                                                                        
  template<typename T>
  Cbasis<T>::~Cbasis() {}


// Assignment Operator          

template<typename T>
Cbasis<T>& Cbasis<T>::operator=(const Cbasis<T>& rhs) {
  if (&rhs == this)
  return *this;

  unsigned new_dim    = rhs.get_dim();
  double new_constant = rhs.get_constant();
  double new_pownc    = rhs.get_pownc();

  elms.resize(new_dim);
    for (unsigned i=0; i<elms.size(); i++) {
    elms[i].resize(2);
    }

    for (unsigned i=0; i<new_dim; i++) {
      for (unsigned j=0; j<2; j++) {
      elms[i][j] = rhs(i, j);
      }
    }
  dim = new_dim;
  constant = new_constant;
  pownc = new_pownc;
  
  return *this;
}

// Conjugate object indices (and switch for operator product evaluation)                              
template<typename T>
Cbasis<T> Cbasis<T>::cjgate(const Cbasis<T>& op) {
  Cbasis result(dim, 0, constant, pownc);
  
  //conjugate indices
  for (unsigned i=0; i<dim; i++) {
  result(i,0) = this->elms[i][1];
  result(i,1) = this->elms[i][0];  
  }

  //assign correct indices   
  for (unsigned i=0; i<op.get_dim(); i++) {
    for (unsigned j=0; j<result.get_dim(); j++) {   
      if( op.elms[i][0]-100 == result(j,1) ) result(j,1) = result(j,1)+100; 
      if( op.elms[i][1]+100 == result(j,0) ) result(j,0) = result(j,0)-100;
    }
  }

  return result;

}


// Access the individual elements             
                             
template<typename T>
  T& Cbasis<T>::operator()(const unsigned& dims, const unsigned& cols) {
  return this->elms[dims][cols];
}

// Access the individual elements (const)

template<typename T>
  const T& Cbasis<T>::operator()(const unsigned& dims, const unsigned& cols) const {
  return this->elms[dims][cols];
}

//Get dimensions
                                                                        
template<typename T>
  unsigned Cbasis<T>::get_dim() const {
  return this->dim;
}

template<typename T>
  double Cbasis<T>::get_constant() const {
  return this->constant;
}

template<typename T>
  double Cbasis<T>::get_pownc() const {
  return this->pownc;
}


// glue together two color strucs                                   
template<typename T>
Cbasis<T>  Cbasis<T>::cadd( const Cbasis<T>& rhs ) {
  Cbasis result(rhs.elms.size() + elms.size(), 0, constant*rhs.get_constant(), pownc + rhs.get_pownc());
  result.elms.clear();
  result.elms.insert(result.elms.end(),rhs.elms.begin(),rhs.elms.end());
  result.elms.insert(result.elms.end(),elms.begin(),elms.end());
  return result;
}


//Print a Cbasis element useful for full Nc Result
template <typename T>
void Cbasis<T>::printBasInfo() 
{
  std::cout << constant << "*Nc^" << pownc << " " << std::endl;
    for (int i=0; i<2; i++) { 
    std::cout << "[ ";
      for(int j=0; j<dim; j++){
      if( elms[j][i] > 0) std::cout << " ";
      if( fabs(elms[j][i]) < 100) std::cout << "  ";
      std::cout << elms[j][i] << " ";
      }
    std::cout << "]";
    std::cout << std::endl;
    }
 std::cout << std::endl;
}


//Fill cbasis indices
template<typename T>
void Cbasis<T>::fill4(const int& i, const int& j, const int& k,const int& l){
    this->elms[0][0]=i;
    this->elms[0][1]=j;
    this->elms[1][0]=k;
    this->elms[1][1]=l;    
}

template<typename T>
void Cbasis<T>::fill6(const int& i, const int& j, const int& k,const int& l,const int& m,const int& n){
    this->elms[0][0]=i;
    this->elms[0][1]=j;
    this->elms[1][0]=k;
    this->elms[1][1]=l;    
    this->elms[2][0]=m;
    this->elms[2][1]=n;    
}


template<typename T>
void Cbasis<T>::fill8(const int& i, const int& j, const int& k, const int& l, 
		      const int& m, const int& n, const int& o,const int& p){
    this->elms[0][0]=i;
    this->elms[0][1]=j;
    this->elms[1][0]=k;
    this->elms[1][1]=l;    
    this->elms[2][0]=m;
    this->elms[2][1]=n;    
    this->elms[3][0]=o;
    this->elms[3][1]=p;    
}

template<typename T>
void Cbasis<T>::app(const double& newCon, const double& npnc){
  this->constant=newCon;
  this->pownc=npnc;
}


// end of member functions

//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------

// evaluates color product <c|T_i T_j|c> for either 
// fundemental or delta representation 

inline std::vector< std::vector< double > > 
colProT( std::vector< std::vector< Cbasis<int> > > c, 
	    std::vector< Cbasis<int> > op)
{
  double hold;
  unsigned test;
  Cbasis<int> take(0,0,0,0); 
  std::vector< std::vector<double> > mat( c.size() );
  
  for(unsigned i = 0; i<c.size(); i++){
    mat[i].resize(c.size());
    for(unsigned j = i; j<c.size(); j++){
      mat[j].resize(c.size());
      for(unsigned k = 0; k<c[i].size(); k++){
        for(unsigned l = 0; l<c[j].size(); l++){
          for(unsigned m = 0; m<op.size();  m++){
	    take = op[m].cadd(c[j][l].cjgate(op[m]));
	    take = c[i][k].cadd(take);
	    ATOOLS::Expression expression;
	    expression.clear();
	    unsigned n = 0;
	    unsigned o = 0;
	    unsigned u = 0;
            while(n < c[i][k].get_dim()){
	      expression.push_back(ATOOLS::Fundamental::New(
	      //take.elms[n][1],fabs(take.elms[n][0]),fabs(take.elms[n][1])));
	      expression.push_back(ATOOLS::Delta::New(take.elms[n][0],take.elms[n][1]));
	    n++;
	    }
	    while(o < op[m].get_dim()){
	      expression.push_back(ATOOLS::Fundamental::New(
	      //take.elms[n+o+u][0],n+fabs(take.elms[n+o][0]),n+fabs(take.elms[n+o][1])));
	      expression.push_back(ATOOLS::Delta::New(take.elms[n+o][0],take.elms[n+o][1]));
	    o++;
	    }
	    while(u < c[j][l].get_dim()){
	      expression.push_back(ATOOLS::Fundamental::New(
	      //take.elms[n+o+u][0],n+o+fabs(take.elms[n+o+u][0]),n+o+fabs(take.elms[n+o+u][1])));
	      expression.push_back(ATOOLS::Delta::New(take.elms[n+o+u][0],take.elms[n+o+u][1]));
	    u++;
	    }
	    expression.Evaluate();

	    //build Nc = 3 matrix
	    hold = take.get_constant()*
	      pow(constants::Nc,take.get_pownc())*real(expression.Result());
	    mat[i][j] += hold; 
            if(i!=j) mat[j][i] += hold; 
	    }
	 }
       }
     }
   }
  return mat;
}


// returns set of C-basis vectors (permuted over a single partition)
inline std::vector< std::vector< Cbasis<int> > > getSinPerms(int numglu,int start)
  {
  std::vector< std::vector< Cbasis<int> > > results;
  double norm  = 1;
  Cbasis<int>  cele(numglu,0,1,0); 
  Cbasis<int>  celeR(numglu,0,1,0);
  
  //Fill standard flow (12...N) + reflection (1N...2)
  if(numglu == 2) {
      for(unsigned j = 0; j<numglu; j++){
         cele.elms[j][0]=start+j+1;
       	 cele.elms[j][1]=-(start+(j+1)%(numglu)+1);
         }
         results.resize(1);
         results[0].push_back(cele);
         return results;
         }
   else for(unsigned j = 0; j<numglu; j++){
         cele.elms[j][0]= start+j+1;
         cele.elms[j][1]= -(start+(j+1)%(numglu)+1);
         celeR.elms[(numglu-1)-j][1]= -(start+j+1);
         celeR.elms[(numglu-1)-j][0]= start+(j+1)%(numglu)+1;
   }
 
    //Generate permutations
   int perms = factorial(numglu-1)/2-1;
   results.resize(perms+1);  
      for(unsigned i = 0; i<=perms; i++){
        results[i].push_back(cele);
        results[i].push_back(celeR);
        xorSwap(&cele.elms[numglu-i%(numglu-2)-1][0], &cele.elms[numglu-i%(numglu-2)-2][0]);
        xorSwap(&cele.elms[numglu-i%(numglu-2)-2][1], &cele.elms[numglu-i%(numglu-2)-3][1]);
        xorSwap(&celeR.elms[i%(numglu-2)+1][0], &celeR.elms[i%(numglu-2)+2][0]);
        xorSwap(&celeR.elms[i%(numglu-2)][1], &celeR.elms[i%(numglu-2)+1][1]); 
   }
  
   return results;
   } 

// Permutes disconnected color objects
inline std::vector< std::vector< Cbasis<int> > >  
  PermT( std::vector< Cbasis<int> > rhsF ,  std::vector< Cbasis<int> >  rhsR){
  std::vector< std::vector< Cbasis<int> > > results;
  int lhsize = rhsF[0].get_dim(); 
  int rhsize = rhsR[0].get_dim();
  int size = lhsize + rhsize;
  int perms = factorial(size)/(factorial(lhsize)*factorial(size-lhsize));
  results.resize(perms);  
  
  for(unsigned k = 0; k<perms; k++){
    //for(unsigned j = 0; j<rhsR.size(); j++){ 
      //for(unsigned i = 0; i<rhsF.size(); i++){
       	 
       results[k].push_back((rhsR[0]).cadd(rhsF[0]));
       //results[k].push_back((rhsR[1]).cadd(rhsF[i]));
       xorSwap(&rhsF[0].elms[1][0], &rhsR[0].elms[k%rhsR[0].get_dim()][0]);
       xorSwap(&rhsF[0].elms[0][1], &rhsR[0].elms[(k+rhsR[0].get_dim()-1)%rhsR[0].get_dim()][1]);  
       //  }
    }
  //}
 

  return results;
}
 

 // Arrange permutations over different color partitions
 inline std::vector< std::vector< Cbasis<int> > > 
 PermPart(std::vector < std::vector< std::vector< Cbasis<int> > > > rhs){
 std::vector< std::vector< Cbasis<int> > > results, hold, take;
 results = rhs[0];
 
 if(rhs.size() == 1 ) return results;

 //Permutations
 for(unsigned i = 1; i < rhs.size(); i++){
   hold.clear();
   for(unsigned l = 0; l < results.size(); l++){
      for(unsigned j = 0; j < rhs[i].size(); j++){
       take = PermT(results[l],rhs[i][j]); 
       hold.insert(hold.end(), take.begin(), take.end());
       }
    }
  results = hold;
  }
 
return results;
}


//returns a complete basis
 inline std::vector< std::vector< Cbasis<int> > > 
  getCbas(std::vector< std::vector<int> > numglu, int nquark, int naquark){
  std::vector< std::vector< Cbasis<int> > > results, test;
  std::vector<std::vector< std::vector< Cbasis<int> > > > iresults;
  iresults.clear();

    for(unsigned g = 0; g < numglu.size(); g++){
      
      //generate single structure permuations
      int count = 0;
      for(unsigned h = 0; h < numglu[g].size(); h++){
      iresults.push_back(getSinPerms(numglu[g][h],count)); 
      count = numglu[g][h];
      }

      //Permute partitions
      test = PermPart(iresults);
      results.insert(results.end(), test.begin(), test.end());
      iresults.clear(); 
    
    }
    return results;
}


// returns minus a matrix
inline void Minus(std::vector< std::vector <double> > &ts){ 
  for(unsigned i = 0; i<ts.size(); i++){
    for(unsigned j = 0; j<ts[i].size(); j++){
    ts[i][j] = -ts[i][j];
    }
  }
}


//Print a monomial
template <typename T>
inline void printMon(std::vector< std::vector <T> > &rhs){
  for(unsigned i = 0; i<rhs.size(); i++){
  std::cout << rhs[i][0] << std::endl;
  std::cout << rhs[i][1] << std::endl;
  }
  std::cout << std::endl;
}




//Print a matrix
template <typename T>
void printMat(std::vector< std::vector <T> > &rhs) {
  std::cout << std::endl;
  for (unsigned i=0; i<rhs.size(); i++) {
    for (unsigned j=0; j<rhs[i].size(); j++) {
    std::cout << rhs[i][j] << ", ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}



#endif
