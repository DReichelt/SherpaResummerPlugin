#ifndef __CBASIS_H
#define __CBASIS_H

#include <vector>
#include <cstdio>
#include <iostream>

const double s_Nc = 3;
const double s_CF = (s_Nc*s_Nc-1.)/(2.*s_Nc);
const double s_CA = s_Nc;
const double s_TR = 1./2.;
const double s_eps = .0001;

using namespace std;

template <typename T> class CBasis
{
private:
    double constant;
    double pownc;

public:
    CBasis(double constant, double pownc);
    CBasis(const CBasis<T>& rhs);
    virtual ~CBasis();

    std::vector< std::vector<T> > Tin;
    std::vector< std::vector<T> > Din;
    std::vector< std::vector<T> > Fin;

    CBasis<T>& operator=(const CBasis<T>& rhs);

    template<typename L>
    friend bool operator==(const CBasis<L>& lhs,const CBasis<L>& rhs);

    CBasis<T>& clear();

    //Glue together two basis elements
    CBasis<T> cadd( const CBasis<T>& rhs );

    //Add color tensors
    CBasis<T>& Tadd(const T& adj,const T& fun,const T& afun);
    CBasis<T>& Dadd(const T& fun,const T& afun);
    CBasis<T>& Fadd(const T& adj1,const T& adj2,const T& adj3);

    //Conjugate basis element
    CBasis<T> cjgate(const CBasis<T>& op);

    //extract perms in COMIX form
    std::vector<int> extract();

    //Append the constants and power of Nc
    void app(const double& newCon, const double& npnc);

    unsigned get_Tdim() const;
    unsigned get_Ddim() const;
    unsigned get_Fdim() const;
    double get_constant() const;
    double get_pownc() const;

    void printBasInfo();


};

typedef std::vector< std::vector< CBasis<int> > > abasis;
typedef std::vector< std::vector< int > > Int_mat;

#endif
