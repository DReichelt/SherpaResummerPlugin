#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream>
//Helpers

void print(int n, int * a) {
    int i ; 
    for (i = 0; i <= n; i++) {
        printf("%d", a[i]); 
    }
    printf("\n"); 
}

inline std::vector< int > cop(int n, int * a) {
  std::vector< int > ires;
  ires.clear();
    for (unsigned i = 0; i <= n; i++) {
    ires.push_back(a[i]);     
    }
    return ires;
}

void integerPartition(int n, int * a, int level, std::vector< std::vector<int> >& rhs){
    int first; 
    int i;
    if (n < 1) return; 
        a[level] = n;
	//print(level, a);
	rhs.push_back(cop(level, a));
	first = (level == 0) ? 2 : a[level-1];
	for(i = first; i <= n / 2; i++){
        a[level] = i; 
        integerPartition(n - i, a, level + 1, rhs);
	}
}


int main(int argc, char ** argv){
    std::vector< std::vector<int> > part;
    int ng = 4;     
    int * a = (int * ) malloc(sizeof(int) * ng); 
    integerPartition (ng, a, 0, part); 

    //Construct large Nc C_basis elements
    


    for(unsigned i = 0; i < part.size(); i++){
      for(unsigned j = 0; j < part[i].size(); j++){
    	std::cout << part[i][j] << " ";
        }
        std::cout << std::endl;
    }

return(0);

}
