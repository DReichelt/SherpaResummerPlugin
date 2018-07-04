#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream>

//Helpers for returning vector of integer partitions
void xorSwap (int* x, int* y) {
    if (x != y) { //ensure that memory locations are different
       *x ^= *y;
       *y ^= *x;
       *x ^= *y;
    }
}

unsigned int factorial(unsigned int n) 
{
    if (n == 0)
       return 1;
    return n * factorial(n - 1);
}

void print(int n, int * a) {
    int i ; 
    for (i = 0; i <= n; i++) {
        printf("%d", a[i]); 
    }
    printf("\n"); 
}


void ppart(std::vector< std::vector<int> > rhs){
for(unsigned i = 0; i < rhs.size(); i++){
      for(unsigned j = 0; j < rhs[i].size(); j++){
    	std::cout << rhs[i][j] << " ";
        }
        std::cout << std::endl;
    }
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
