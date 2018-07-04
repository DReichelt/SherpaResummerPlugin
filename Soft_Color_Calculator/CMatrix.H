#ifndef __CMatrix_H
#define __CMatrix_H


class CMatrix
{
private:

public:
    Complex* mat;
    CMatrix();
    ~CMatrix();
    CMatrix& operator=(const CMatrix& rhs);
    void printMat(int dim);
    void H_conjugate(int dim);
    void Minus(int dim);
    void printFile(std::ofstream& file, int dim);
    void colProT( abasis c , std::vector< CBasis<int> > Coperator, double v_Nc = 3 );
    void cleanup();
};
void print_cl_info(std::string fn);
CBasis<int> cjgate(const CBasis<int>& cconj, const CBasis<int>& op);
inline std::vector< CBasis<int> > Ciden(const int& in, const int& ln);
inline std::vector< CBasis<int> > MC(const int& in, const int& jn, const int& kn,const int& ln);
inline std::vector< CBasis<int> > MPC(const int& io, const int& jo, const int& it, const int& jt);
inline std::vector< CBasis<int> >
PC(const int& A, const int& B, const int& D, const int& E);
inline std::vector< CBasis<int> > PCT(const int& io, const int& it, const int& ith, const int& ifo);
void Texpre(CBasis<int> ct, ATOOLS::Expression &ex, int &pl);
void normalize_basis(abasis &unc, double v_Nc = 3);
void compute_metric(CMatrix& result, abasis c);
void compute_color_product(abasis NBasis, std::vector< CMatrix >& T, int nq, int ng, int numprods, double Nc);
void check_color_conservation(const std::vector< CMatrix >& T ,const CMatrix Cmetric,const int ng,const int nq,const int bdim);

#endif