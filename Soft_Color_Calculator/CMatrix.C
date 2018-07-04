#ifndef __CMatrix_C
#define __CMatrix_C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstdio>
#include <vector>
#include <iostream>
#include "CBasis.H"
#include "CBasis.C"
#include "CMatrix.H"

CMatrix::CMatrix() {}
CMatrix::~CMatrix() {}


CMatrix& CMatrix::operator=(const CMatrix& rhs)
{
    if (&rhs == this)
        return *this;
    return *this;
}


//Print a matrix
void CMatrix::printMat(int dim)
{
    std::cout << std::endl;
    for (unsigned i=0; i<dim; i++)
    {
        std::cout << "{";
        for (unsigned j=0; j<dim; j++)
        {
            std::cout << std::setprecision(15) << this->mat[i*dim+j].real();
            if(j<dim-1)
            {
                cout <<"," ;
            }
        }
        std::cout <<((i+1==dim)?"}\n":"},\n");
    }
    std::cout <<"}]\n";
    std::cout << "===================================" << std::endl;
    std::cout << std::endl;
}


//Conjugate matrices
void CMatrix::H_conjugate(int dim)
{
    complex<double> hold;
    for (size_t i=0; i<dim; i++)
    {
        for (size_t j=i; j<dim; j++)
        {
            hold = this->mat[i*dim+j];
            this->mat[i*dim+j] = conj(this->mat[j*dim+i]);
            this->mat[j*dim+i] = conj(hold);
        }
    }
}

// returns minus a matrix
void CMatrix::Minus(int dim)
{
    for(unsigned i = 0; i<dim; i++)
    {
        for(unsigned j = 0; j<dim; j++)
        {
            this->mat[i*dim+j] = -this->mat[i*dim+j];
        }
    }
}

//Print to file
void CMatrix::printFile(std::ofstream& file, int dim)
{
    file << std::endl;
    for (unsigned i=0; i<dim; i++)
    {
        for (unsigned j=0; j<dim; j++)
        {
            file << std::setprecision(15) << this->mat[i*dim+j].real() << " ";
        }
        //file << "\n";
    }
    //file <<"\n";
}


void CMatrix::colProT( abasis c , std::vector< CBasis<int> > Coperator, double v_Nc )
{
    double hold;
    int m_dim = c.size();
    CBasis<int> take(1,0);
    this->mat = new Complex[m_dim*m_dim];

    for(unsigned i = 0; i<m_dim; i++)
    {
        for(unsigned j = i; j<m_dim; j++)
        {
            for(unsigned k = 0; k<c[i].size(); k++)
            {
                for(unsigned l = 0; l<c[j].size(); l++)
                {
                    for(unsigned m = 0; m<Coperator.size(); m++)
                    {
                        ATOOLS::Expression expression;
                        expression.clear();
                        int place = 0;

                        Texpre(c[i][k], expression, place);
                        Texpre(Coperator[m], expression, place);
                        Texpre(cjgate(c[j][l],Coperator[m]), expression, place);

                        //output color string to be evaluated

                        expression.Evaluate();


                        //build Nc = 3 matrix
                        hold = real(expression.Result()) + imag(expression.Result());
                        this->mat[i*m_dim+j] += hold;
                        if( fabs( this->mat[i*m_dim+j].real() ) < s_eps )
                        {
                            hold = 0;
                            this->mat[i*m_dim+j]=0;
                            this->mat[j*m_dim+i]=0;
                        }
                        if(i!=j) this->mat[j*m_dim+i] += hold;
                    }
                }
            }
        }
    }
}

void CMatrix::cleanup()
{
    delete[] this->mat;
    return;
}


// Complex conjugates the indices and
// decides which indices the exchange operator
// is speaking to...in which case it adjusts the
// indices accordingly by adding 100.

void print_cl_info(std::string fn)
{
    std::cout << std::endl;
    std::cout << "Process: " << fn << std::endl;
    std::cout << std::endl;
    std::cout << "===================================" << std::endl;
    std::cout << "Color basis manipulator" << std::endl;
    std::cout << "===================================" << std::endl;
    std::cout << std::endl;
    return;
}

CBasis<int> cjgate(const CBasis<int>& cconj, const CBasis<int>& op)
{
    CBasis<int> results(1,0);
    results = cconj;

    //Conjugate the fundemental indices
    for( unsigned i = 0; i<cconj.get_Ddim(); i++)
    {
        results.Din[i][0] = fabs(cconj.Din[i][1]);
        results.Din[i][1] = fabs(cconj.Din[i][0]);
    }

    //Conjugate the T fund indices
    for( unsigned i = 0; i<cconj.get_Tdim(); i++)
    {
        results.Tin[i][1] = fabs(cconj.Tin[i][2]);
        results.Tin[i][2] = fabs(cconj.Tin[i][1]);
    }

    //assign fun-antifundamental indices T <-> F/T
    for (unsigned i=0; i<op.get_Ddim(); i++)
    {
        for (unsigned j=0; j<cconj.get_Tdim(); j++)
        {
            if( op.Din[i][1]-100 == results.Tin[j][1] ) results.Tin[j][1] = results.Tin[j][1]+100;
            if( op.Din[i][0]-100 == results.Tin[j][2] ) results.Tin[j][2] = results.Tin[j][2]+100;
        }
        for (unsigned j=0; j<cconj.get_Ddim(); j++)
        {
            if( op.Din[i][0]-100 == results.Din[j][1] ) results.Din[j][1] = results.Din[j][1]+100;
            if( op.Din[i][1]-100 == results.Din[j][0] ) results.Din[j][0] = results.Din[j][0]+100;
        }
    }

    //assign fun-antifundamental indices T <-> F/T
    for (unsigned i=0; i<op.get_Tdim(); i++)
    {
        for (unsigned j=0; j<cconj.get_Tdim(); j++)
        {
            if( op.Tin[i][2]-100 == results.Tin[j][1] && results.Tin[j][1]>50) results.Tin[j][1] = results.Tin[j][1]+100;
            if( op.Tin[i][1]-100 == results.Tin[j][2] && results.Tin[j][2]>50) results.Tin[j][2] = results.Tin[j][2]+100;
        }
        for (unsigned j=0; j<cconj.get_Ddim(); j++)
        {
            if( op.Tin[i][1]-100 == results.Din[j][1] ) results.Din[j][1] = results.Din[j][1]+100;
            if( op.Tin[i][2]-100 == results.Din[j][0] ) results.Din[j][0] = results.Din[j][0]+100;
        }
    }

    //assign correct adjoint indices T <-> F/T
    for (unsigned i=0; i<op.get_Tdim(); i++)
    {
        for (unsigned j=0; j<cconj.get_Tdim(); j++)
        {
            if( op.Tin[i][0]-100 == results.Tin[j][0] ) results.Tin[j][0] = results.Tin[j][0]+100;
        }
        for (unsigned k=0; k<cconj.get_Fdim(); k++)
        {
            for(unsigned adj = 0; adj < 3; adj++)
            {
                if( op.Tin[i][0]-100 == results.Fin[k][adj] ) results.Fin[k][adj] = results.Fin[k][adj]+100;
            }
        }
    }

    //avoid double counting internal f-indices
    for (unsigned k=0; k<cconj.get_Fdim(); k++)
    {
        for(unsigned adj = 0; adj < 3; adj++)
        {
            if( results.Fin[k][adj] > 500 ) results.Fin[k][adj] = results.Fin[k][adj]+100;
        }
    }

    //assign indices F <-> F
    for (unsigned i=0; i<op.get_Fdim(); i++)
    {
        for (unsigned j=0; j<cconj.get_Fdim(); j++)
        {
            for(unsigned adj = 0; adj < 3; adj++)
            {
                if( op.Fin[i][adj]-100 == results.Fin[j][0]) results.Fin[j][0] = results.Fin[j][0]+100;
                if( op.Fin[i][adj]-100 == results.Fin[j][1]) results.Fin[j][1] = results.Fin[j][1]+100;
                if( op.Fin[i][adj]-100 == results.Fin[j][2]) results.Fin[j][2] = results.Fin[j][2]+100;
            }
        }
        for (unsigned j=0; j<cconj.get_Tdim(); j++)
        {
            for(unsigned adj = 0; adj < 3; adj++)
            {
                if( op.Fin[i][adj]-100 == results.Tin[j][0]) results.Tin[j][0] = results.Tin[j][0]+100;
            }
        }
    }

    return results;
}

//identity operator
inline std::vector< CBasis<int> > Ciden(const int& in, const int& ln)
{
    std::vector< CBasis<int> > results;
    CBasis<int> p1(1,0);
    results.push_back(p1);
    return results;
}


//Exchange operator between q and q
inline std::vector< CBasis<int> > MC(const int& in, const int& jn, const int& kn,const int& ln)
{
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

//Exchange operator between g and q
inline std::vector< CBasis<int> > MPC(const int& io, const int& jo, const int& it, const int& jt)
{
    std::vector< CBasis<int> > results;

    CBasis<int> p1(1,0);
    p1.Tadd(jo,it,25);
    p1.Tadd(io,25,jt);

    CBasis<int> p2(-1,0);
    p2.Tadd(io,it,25);
    p2.Tadd(jo,25,jt);

    results.push_back(p1);
    results.push_back(p2);

    return results;
}

//Exchange operator between g and g in terms of adjoint matrices (faster)
//COMIX knows about the U(1) decoupling terms
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

//Exchange operator between g and g in terms of Fundamental matrices (slower)
inline std::vector< CBasis<int> > PCT(const int& io, const int& it, const int& ith, const int& ifo)
{
    std::vector< CBasis<int> > results;
    CBasis<int> p1(-2.,0);
    p1.Tadd(io,io,it);
    p1.Tadd(it,it,ith);
    p1.Tadd(ith,ith,ifo);
    p1.Tadd(ifo,ifo,io);

    CBasis<int> p2(-2.,0);
    p2.Tadd(io,io,ifo);
    p2.Tadd(it,it,io);
    p2.Tadd(ith,ith,it);
    p2.Tadd(ifo,ifo,ith);

    CBasis<int> p3(2.,0);
    p3.Tadd(io,io,ith);
    p3.Tadd(it,it,io);
    p3.Tadd(ith,ith,ifo);
    p3.Tadd(ifo,ifo,it);

    CBasis<int> p4(2.,0);
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


//write to expression T's computed with "Color" library
//numbering of indices assures no inconsitencies
void Texpre(CBasis<int> ct, ATOOLS::Expression &ex, int &pl)
{
    int pl2 = 0;
    int tA, tf, taf;
    ex.push_back(ATOOLS::CNumber::New(ct.get_constant()*pow(s_Nc,ct.get_pownc())));

    for(unsigned i = 0; i < ct.get_Tdim(); i++)
    {
        tA = ct.Tin[i][0];
        //Decides whether quark index is internal or not
        //external indices are 50-60 or 150-160 depending on whether they
        //are on left or right (conjugate) side of <c| P |c>
        tf  = ((fabs(ct.Tin[i][1]-50) > 10 && fabs(ct.Tin[i][1]-150) > 10) ? 1000+pl : 0) + fabs(ct.Tin[i][1]);
        taf = ((fabs(ct.Tin[i][2]-50) > 10 && fabs(ct.Tin[i][2]-150) > 10) ? 1000+pl : 0) + fabs(ct.Tin[i][2]);
        ex.push_back(ATOOLS::Fundamental::New(tA,tf,taf));
        pl2++;
    }

    for(unsigned i = 0; i < ct.get_Ddim(); i++)
    {
        ex.push_back(ATOOLS::Delta::New(fabs(ct.Din[i][0]),fabs(ct.Din[i][1])));
    }

    for(unsigned i = 0; i < ct.get_Fdim(); i++)
    {
        ex.push_back(ATOOLS::Adjoint::New(ct.Fin[i][0],ct.Fin[i][1],ct.Fin[i][2]));
    }
    //makes sure that internal indices are not double counted
    pl = 2*(pl+pl2)+10;
}


// Normalize the basis
void normalize_basis(abasis &unc, double v_Nc)
{

    double hold,normfac;

    for(unsigned i = 0; i<unc.size(); i++)
    {
        normfac=0;
        for(unsigned k = 0; k<unc[i].size(); k++)
        {
            for(unsigned l = 0; l<unc[i].size(); l++)
            {
                ATOOLS::Expression expression;
                expression.clear();
                //expression.setNC(v_Nc);
                int place = 0;

                Texpre( unc[i][k],expression , place );
                Texpre( cjgate(unc[i][l], Ciden(0,0)[0]) , expression , place );
                expression.Evaluate();

                //append the i-th element
                hold = real(expression.Result());
                normfac = normfac + hold;
            }
        }
        for(unsigned m = 0; m < unc[i].size(); m++)
        {
            unc[i][m].app(unc[i][m].get_constant()*(1/sqrt(normfac)),unc[i][m].get_pownc());
        }
    }
    return;
}

void compute_metric(CMatrix& result, abasis c)
{
    result.colProT(c,Ciden(0,0));
    return;
}


//compute the color products
void compute_color_product(abasis NBasis, std::vector< CMatrix >& T, int nq, int ng, int numprods, double Nc)
{
    int numlegs = ng + 2*nq;
    int naq = nq;
    T.resize(numprods);
    int coin = 0;
    int i_f, j_f;
    int ik, jk;
    int bdim = NBasis.size();
    for(unsigned i = 1; i<=numlegs; i++)
    {
        for(unsigned j = i+1; j<=numlegs; j++)
        {
            i_f = 0;
            j_f = 0;
            ik = i;
            jk = j;
            //Decide exchange type
            //===============================
            //decide if i/j is q or g
            if(i <= nq+naq)
            {
                i_f  = (i <= nq) ? -1 : 1;
            }
            if(j <= nq+naq)
            {
                j_f  = (j <= nq) ? -1 : 1;
            }
            //===============================

            //qq MC(a,b,c,d)
            //===============================
            if(fabs(i_f) == 1 && fabs(j_f) == 1)
            {
                T[coin].colProT( NBasis, MC((i_f==1 ?  50+i  : 100+50+i), (i_f==1 ? 100+50+i : 50+i ),
                                            (j_f==1 ?  50+j  : 100+50+j), (j_f==1 ? 100+50+j : 50+j )   ), Nc );
            }
            //===============================

            //qg MPC(a,b,c,d)
            //===============================
            if(fabs(i_f) == 1 && j_f == 0)
            {
                T[coin].colProT( NBasis, MPC(j - (nq+naq),100+(j - (nq+naq)),
                                             (i_f==1 ? 50+i : 100+50+i  ) , (i_f==1 ? 100+50+i : 50+i)  ), Nc ) ;
            }
            //===============================

            //gg PCT(a,b,c,d)
            //===============================
            if(j_f == 0 && i_f == 0)
            {
                T[coin].colProT( NBasis, PC(100+(i - (nq+naq)),(i - (nq+naq)),(j - (nq+naq)),100+(j - (nq+naq) ) ), Nc );
            }
            //===============================

            // Very Important!!! minus signs for incoming quark, can
            // assign custom signs to external legs through eg:
            // ----------------------------------------------------
            for(int k = 1; k < nq+1; k++)
            {
                if( i == k || j == k ) T[coin].Minus(bdim);
            }
            // ----------------------------------------------------
            coin++;

            i = ik;
            j = jk;
        }
    }
}


//--------------------------------------------
//Check color conservation
void check_color_conservation(const std::vector< CMatrix >& T ,const CMatrix Cmetric,const int ng,const int nq,const int bdim)
{

    if(T.size() == 0) return;

    int naq = nq;
    int numlegs = ng + 2*nq;
    int coin = 0;

    double check;
    double Casmir;
    int k_t,i_t, j_t, loop,count;
    for(unsigned k = 0; k<numlegs; k++)
    {
        std::cout << "Should be T" << k+1 << ".J = 0" << std::endl;
        k_t = (k < 2 ? k : 2*k-1);
        Casmir = k < (nq+naq) ? 4./3. : 3 ;
        if(Casmir == 3) std::cout << "is a " << "g" << std::endl;
        if(Casmir == 4./3.) std::cout << "is a " << "q" << std::endl;
        for(unsigned i = 0; i<bdim; i++)
        {
            for(unsigned j = 0; j<bdim; j++)
            {
                check = real(Casmir*Cmetric.mat[ i*bdim + j]);
                i_t=0;
                j_t=0;
                loop=0;
                count = 0;
                for(unsigned l = 0; l<T.size(); l++)
                {
                    i_t = loop;
                    j_t = count+loop+1;
                    count++;
                    if(j_t == numlegs-1)
                    {
                        loop++;
                        count = 0;
                    }
                    if(i_t == k || j_t == k)
                    {
                        check += real(T[l].mat[i*bdim+j]);
                    }
                }
                std::cout << (fabs(check) > .0001 ? check : 0) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    return;
}


#endif

