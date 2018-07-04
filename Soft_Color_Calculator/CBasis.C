#ifndef __CBASIS_C
#define __CBASIS_C

#include "CBasis.H"

// Parameter Constructor
template<typename T>
CBasis<T>::CBasis(double _constant, double _pownc)
{
    constant = _constant;
    pownc = _pownc;
}

// Copy Constructor
template<typename T>
CBasis<T>::CBasis(const CBasis<T>& rhs)
{
    constant  = rhs.get_constant();
    pownc  = rhs.get_pownc();
    Tin = rhs.Tin;
    Din = rhs.Din;
    Fin = rhs.Fin;
}

// (Virtual) Destructor
template<typename T>
CBasis<T>::~CBasis() {}

// Assignment Operator
template<typename T>
CBasis<T>& CBasis<T>::operator=(const CBasis<T>& rhs)
{
    if (&rhs == this)
        return *this;

    constant = rhs.get_constant();
    pownc    = rhs.get_pownc();

    unsigned new_T_dim = rhs.get_Tdim();
    unsigned new_D_dim = rhs.get_Ddim();
    unsigned new_F_dim = rhs.get_Fdim();

    Tin.resize(new_T_dim);
    for (unsigned i=0; i<Tin.size(); i++)
    {
        Tin[i].resize(3);
        for (unsigned j=0; j<Tin[i].size(); j++)
        {
            Tin[i][j] = rhs.Tin[i][j];
        }
    }

    Din.resize(new_D_dim);
    for (unsigned i=0; i<Din.size(); i++)
    {
        Din[i].resize(2);
        for (unsigned j=0; j<Din[i].size(); j++)
        {
            Din[i][j] = rhs.Din[i][j];
        }
    }

    Fin.resize(new_F_dim);
    for (unsigned i=0; i<Fin.size(); i++)
    {
        Fin[i].resize(3);
        for (unsigned j=0; j<Fin[i].size(); j++)
        {
            Fin[i][j] = rhs.Fin[i][j];
        }
    }

    return *this;
}


// comparison operator
template<typename L>
bool operator==(const CBasis<L>& lhs,const CBasis<L>& rhs)
{
    if((lhs.get_Tdim() != rhs.get_Tdim()) || (lhs.get_Ddim() != rhs.get_Ddim()) || (lhs.get_Fdim() != rhs.get_Fdim())) return 0;

    for(size_t i = 0; i < lhs.get_Tdim(); i++)
    {
        if((lhs.Tin[i][0] != rhs.Tin[i][0]) || (lhs.Tin[i][1] != rhs.Tin[i][1]) || (lhs.Tin[i][2] != rhs.Tin[i][2])) return 0;
    }
    for(size_t i = 0; i < lhs.get_Ddim(); i++)
    {
        if( (lhs.Din[i][0] != rhs.Din[i][0]) || (lhs.Din[i][1] != rhs.Din[i][1])) return 0;;
    }
    for(size_t i = 0; i < lhs.get_Fdim(); i++)
    {
        if( (lhs.Fin[i][0] != rhs.Fin[i][0]) || (lhs.Fin[i][1] != rhs.Fin[i][1]) || (lhs.Fin[i][2] != rhs.Fin[i][2])) return 0;
    }
    return 1;

}

template<typename T>
std::vector<int> CBasis<T>::extract()
{
    std::vector< int > intperms;
    return intperms;
}

template<typename T>
unsigned CBasis<T>::get_Tdim() const
{
    return this->Tin.size();
}
template<typename T>
unsigned CBasis<T>::get_Ddim() const
{
    return this->Din.size();
}
template<typename T>
unsigned CBasis<T>::get_Fdim() const
{
    return this->Fin.size();
}

template<typename T>
double CBasis<T>::get_constant() const
{
    return this->constant;
}

template<typename T>
double CBasis<T>::get_pownc() const
{
    return this->pownc;
}

template<typename T> CBasis<T>& CBasis<T>::clear()
{
    Tin.clear();
    Din.clear();
    Fin.clear();
    constant = 1.;
    pownc = 0;
    return *this;
}

template<typename T>
CBasis<T>& CBasis<T>::Tadd(const T& adj,const T& fun,const T& afun)
{
    std::vector<T> tares(3);
    tares[0] = adj;
    tares[1] = fun;
    tares[2] = afun;
    Tin.push_back(tares);
    return *this;
}

template<typename T>
CBasis<T>& CBasis<T>::Dadd(const T& fun,const T& afun)
{
    std::vector<T> tares(2);
    tares[0] = fun;
    tares[1] = afun;
    Din.push_back(tares);
    return *this;
}

template<typename T>
CBasis<T>& CBasis<T>::Fadd(const T& adj1,const T& adj2,const T& adj3)
{
    std::vector<T> tares(3);
    tares[0] = adj1;
    tares[1] = adj2;
    tares[2] = adj3;
    Fin.push_back(tares);
    return *this;
}

// glue together two color strucs
template<typename T>
CBasis<T> CBasis<T>::cadd( const CBasis<T>& rhs )
{
    CBasis result(constant*rhs.get_constant(), pownc + rhs.get_pownc());
    result.Tin.insert(result.Tin.end(),Tin.begin(),Tin.end());
    result.Tin.insert(result.Tin.end(),rhs.Tin.begin(),rhs.Tin.end());
    result.Din.insert(result.Din.end(),Din.begin(),Din.end());
    result.Din.insert(result.Din.end(),rhs.Din.begin(),rhs.Din.end());
    result.Fin.insert(result.Fin.end(),Fin.begin(),Fin.end());
    result.Fin.insert(result.Fin.end(),rhs.Fin.begin(),rhs.Fin.end());
    return result;
}


//Print a Cbasis element
template <typename T>
void CBasis<T>::printBasInfo()
{
    std::cout << constant << "*Nc^" << pownc << " " << std::endl;
    if(Tin.size() > 0)
    {
        cout << "T.fund" << endl;
        for (int i=0; i<3; i++)
        {
            std::cout << "[ ";
            for(int j=0; j<Tin.size(); j++)
            {
                if( Tin[j][i] > 0) std::cout << " ";
                std::cout << Tin[j][i] << " ";
            }
            std::cout << "]";
            std::cout << std::endl;
        }
    }

    if(Din.size() > 0)
    {
        cout << "D.fund" << endl;
        for (int i=0; i<2; i++)
        {
            std::cout << "[ ";
            for(int j=0; j<Din.size(); j++)
            {
                if( Din[j][i] > 0) std::cout << " ";
                std::cout << Din[j][i] << " ";
            }
            std::cout << "]";
            std::cout << std::endl;
        }
    }

    if(Fin.size() > 0)
    {
        cout << "F.Adjoint" << endl;
        for (int i=0; i<3; i++)
        {
            std::cout << "[ ";
            for(int j=0; j<Fin.size(); j++)
            {
                if( Fin[j][i] > 0) std::cout << " ";
                std::cout << Fin[j][i] << " ";
            }
            std::cout << "]";
            std::cout << std::endl;
        }
    }
    std::cout << std::endl;
}


template<typename T>
void CBasis<T>::app(const double& newCon, const double& npnc)
{
    this->constant=newCon;
    this->pownc=npnc;
}

#endif


