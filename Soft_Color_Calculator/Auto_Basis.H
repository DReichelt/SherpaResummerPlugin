#ifndef __Auto_Basis_H
#define __Auto_Basis_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstdio>
#include <vector>
#include <iostream>
#include "CBasis.H"

template <typename T> class Auto_Basis
{
private:
public:

    std::vector < CBasis<int> > ToFundamental(std::string rhs, int nquarks, std::string qs,  std::string aqs);
    CBasis<int> ToAdjoint(std::string rhs);
    inline CBasis<int> AdjointDelta(int a, int b);
};
abasis auto_Cbasis( int number_of_gluons, int number_of_quarks_pairs, const bool &adjoint=false);
void printBasis(abasis rhs);
inline std::vector< CBasis <int> > combineA(std::vector< CBasis <int> > rhs1 , std::vector< CBasis <int> > rhs2 );
inline abasis fold(abasis rhs1, abasis rhs2);
inline CBasis<int>   AdjointDelta(int a, int b);
inline void  get_partition_symmetry_factor(std::vector< int > partition, size_t &sym_factor, size_t &sym_s);
inline std::vector< std::string > arrange_groups( std::vector<int> parts, std::string labels, bool is_quarks);
abasis TraceConGluons( std::string rhs, int num_quarks, std::string ql);
abasis gluon_groups( int ngluon, std::string gluon_labels, const bool &multip=false );
abasis quark_groups(int nquark, std::string quark_labels);
abasis connected_qqng_groups(int ngluon,int quark_pairs,
                             int con_qp, int con_gl ,std::string gluon_labels,std::string quark_labels);
abasis qqng_groups(int ngluon, int quark_pairs,std::string gluon_labels, std::string quark_labels, const bool &adjoint=false);
void printBasis(abasis rhs);

#endif
