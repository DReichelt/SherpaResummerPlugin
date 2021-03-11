#ifndef MATRIX_H
#define MATRIX_H
/** \file */
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/Random.H"
#include <valarray>
#include <cassert>
#include <complex>

#include "Math/matexp.hpp"
#include "Math/asa007.hpp"
#include "Math/linpack_d.hpp"


/// Basic matrix class.
/**
    Class implementing basic functionality. If you need more, use some
    linear algebra library.
*/

namespace RESUM {

  template <typename T>
  T _conjugate(const T& z) { return z; }

  template <typename X>
  std::complex<X> _conjugate(const std::complex<X>& z) { return std::conj(z); }
  
  template <typename T>
  class Matrix {

  protected:
    /// Valarray with underlying data.
    std::valarray<T> m_data;

    /// Number of columns.
    size_t m_dim;
  
  public:
    /// \name Standard constructors.
    // \{
    /// Default constructor.
    Matrix() = default;
    
    /// Copy constructor.
    Matrix<T>(const Matrix<T>& other) = default;
    
    /// Copy constructor for matrix with arbitrary type (that can be cast to T).
    template <typename OTHER>
    Matrix<T>(const Matrix<OTHER>& other) {
      m_dim = other.numCols();
      m_data = std::valarray<T>(other.numElements());
      for(size_t i=0; i<other.numElements(); i++) m_data[i] = T(other.data()[i]);
    }
    
    /// Move constructor.
    Matrix<T>(Matrix<T>&& other) : Matrix<T>() {
      other.swap(*this);
    }
    
    /// Constructor directly using underlying \p data and number of cols \p c.
    /**
       \param data Valarray with underlying data.
       \param c Number of columns to be extracted.
       
       The number of elements in \p data has to be an integer multiple 
       (the number of rows) of \p c.
    */
    Matrix(const std::valarray<T> data, size_t c) : m_data(data), m_dim(c) {
      assert(data.size()%c == 0);
    }
    
    /// Construct empty matrix with \p r rows and \p c cols.
    /**
       \param r number of rows
       \param c number of columns
       \param def default entry
    */
    Matrix(size_t r, size_t c, T def=T()) : m_data(def,r*c), m_dim(c) {}
    
    /// Construct matrix from nested vectors, determine number of rows and cols from vector.
    /**
       \param mat Nested vector containing underlying data.
       
       Every entry in \p mat will be interprereted as a row of the matrix. 
       The first row determines the number of columns. If one of the later rows is longer the additional elements
       are ignored. If one of the later rows is shorter, an error will be thrown.
    */
    Matrix(const std::vector<std::vector<T> >& mat) {
      _init(mat);
    }
    
    /// Construct matrix from nested vectors, but fill with zeros to match \p dim rows and cols.
    /**
       \param mat Nested vector containing underlying data.
       \param dim Dimension of matrix to be constructed.
       \param def Value to be appended to each row to contain \p dim number of elements
       
       Every entry in \p mat will be interprereted as a row of the matrix. 
       If the rows have different sizes an error will be thrown.
    */
    Matrix(const std::vector<std::vector<T> >& mat, size_t rows,
           size_t cols, T def=T()) {
      _init(mat,rows,cols,def);
    }
    
    
    /// Construct matrix from pointer to nested vectors, determine number of rows and cols from vector.
    /**
       \param mat Pointer to nested vector containing underlying data.
       
       Every entry in \p mat will be interprereted as a row of the matrix. 
       If the rows have different sizes an error will be thrown.
       If \p mat is a nullptr, an empty matrix will be constructed.
    */
    Matrix(const std::vector<std::vector<T> >* mat) {
      if(mat) {
        _init(*mat);
      }
      else {
        m_data = std::valarray<T>();
        m_dim = 0;
      }
    }
    
    /// Construct matrix from pointer to nested vectors, but fill with zeros to match \p dim rows and cols.
    /**
       \param mat Pointer to nested vector containing underlying data.
       \param dim Dimension of matrix to be constructed.
       \param def Value to be appended to each row to contain \p dim number of elements
       
       Every entry in \p mat will be interprereted as a row of the matrix. 
       If the rows have different sizes an error will be thrown.
       If \p mat is a nullptr, a square matrix of with \p dim rows and columns, 
       and every element equal to \p def will be constructed.
    */
    Matrix(const std::vector<std::vector<T> >* mat, size_t rows,
           size_t cols, T def=T()) {
      if(mat) {
        _init(*mat, rows, cols, def);
      }
      else {
        m_data = std::valarray<T>(def,rows*cols);
        m_dim = cols;
      }
    }
    
    // \}
    
    /// Build diagonal square matrix of size \p dim.
    /**
       \param diag Value of diagonal elements.
       \param dim Number of rows and columns of the matrix.
       \param offdiag Value of the offdiagonal elements.
    */
    static Matrix<T> diagonal(const T& diag, size_t dim,
                              const T& offdiag=T()) {
      Matrix<T> ret(dim,dim,offdiag);
      for(size_t i=0; i<dim;i++) ret(i,i) = diag;
      return ret;
    }
    
    /// Build random matrix with \p r rows and \p c cols.
    /**
       \param r number of rows
       \param c number of columns
    */
    static Matrix<T> random(size_t r, size_t c) {
      Matrix<T> ret(r,c);
      for(T& d: ret.m_data) d = ATOOLS::ran->Get();
      return ret;
    }
    
    /// \name Assignments.
    // \{
    /// Swap algorithm.
    void swap (Matrix<T>& other) noexcept {
      std::swap(m_data, other.m_data);
      std::swap(m_dim, other.m_dim);
    }
    
    /// copy and move assignment
    Matrix<T>& operator=(Matrix<T> other) {
      other.swap(*this);
    return *this;
    }
    // \}
    
    /// \name Basic mathematical functionality.
    // \{

    /// Matrix addition.
    Matrix<T>& operator+=(const Matrix<T>& rhs) {
      assert(sameDim(rhs));
      m_data += rhs.m_data;
      return *this;
    }
    
    /// Matrix subtraction.
    Matrix<T>& operator-=(const Matrix<T>& rhs) {
      assert(sameDim(rhs));
      m_data -= rhs.m_data;
      return *this;
    }
    
    /// Matrix multiplication.
    Matrix<T>& operator*=(const Matrix<T>& other) {
      assert(numCols() == other.numRows());
      const size_t rows = numRows();
      const size_t cols = numCols();
      const size_t ocols = other.numCols();
      Matrix<T> tmp(rows,ocols,T(0));
      transposeInPlace();
      this->swap(tmp);
      for(size_t k=0; k<cols; k++) {
        for(size_t i=0; i<rows; i++) {
          for(size_t j=0; j<ocols; j++) {
            element(i,j) += tmp(k,i)*other(k,j);
          }
        }
      }
      return *this;
    }
    
    /// Multiplication with scalar.
    Matrix<T>& operator*=(T b) {
      m_data *= b;
      return *this;
    }
    
    /// Division by scalar.
    Matrix<T>& operator/=(T b) {
      m_data /= b;
      return *this;
    }
    
    /// Calculate the trace of the matrix.
    /**
       Throws an error if matrix is not square.
    */
    T trace() const {
      if(!isSquare()) THROW(fatal_error,"Trace of non-square matrix is undefined!");
      return m_data[std::slice(0, m_dim, m_dim+1)].sum();
    }
    
    /// Transpose the matrix in-place.
    /**
       This is really in-place only for square matrices.
    */
    Matrix<T> transposeInPlace() {
      if(isSquare()) {
        for(size_t r=0; r<numRows(); r++) {
          for(size_t c=r+1; c<numCols(); c++) {
            std::swap(element(r,c), element(c,r));
          }
        }
      }
      else {
        Matrix<T> tmp = *this;
        m_dim = numRows();
        for(size_t r=0; r<numRows(); r++) {
          for(size_t c=0; c<numCols(); c++) {
            element(r,c) = tmp(c,r);
          }
        }
      }
      return *this;
    }
    
    /// Return the transpose of the matrix.
    Matrix<T> transpose() const {
      return Transpose(*this);
    }
    
    /// Calculate the conjugate of the the matrix in-place.
    Matrix<T> conjugateInPlace() {
      if(std::is_same<T,std::complex<double>>::value) {
        for(T& d: m_data) d = _conjugate(d);
      }
      // m_data = m_data.apply(std::conj);
      return *this;
    }
    
    /// Return the hermitian transpose of the matrix.
    Matrix<T> conjugate() const {
      return Conjugate(*this);
    }
    
  
    /// Calculate the hermitian transpose of the the matrix in-place.
    /**
       This is really in-place only for square matrices.
    */
    Matrix<T> hermitianTransposeInPlace() {
      transposeInPlace();
      conjugateInPlace();
      return *this;
    }
    
    /// Return the hermitian transpose of the matrix.
    Matrix<T> hermitianTranspose() const {
      return HermitianTranspose(*this);
    }
    
    /// Calculate matrix exponential of this matrix in place;
    Matrix<T> exponentiateInPlace() {
      T* result = c8mat_expm1(numCols(), &m_data[0]);
      m_data = {result, numElements()};
      delete [] result;
      return *this;
    }

    /// Calculate matrix exponential of this matrix;
    Matrix<T> exponentiate() {
      return exp(*this);
    }
    
    /// Calculate the inverse of the symmetric matrix with lower triangular matrix equal to this one in place.
    Matrix<T> inverseLowerTriangularInPlace() {
      Matrix<T> old(*this);
      _toLowTriang();
      int nullity;
      int ifault;
      std::vector<double> W(m_dim);
      std::valarray<T> new_data(m_dim*(m_dim+1)/2);
      syminv(&m_data[0], m_dim, &new_data[0], &W[0],
             &nullity, &ifault);
      // ifault = 2;
      switch(ifault) {
      case 0:
        msg_Debugging()<<"Rank of deficiency: "<<nullity<<"\n";
        _fromLowTriang(std::move(new_data));
        break;
      case 1: msg_Error()<<"N < 1.\n"; break;
      case 2:
        {
          msg_Debugging()<<"A is not positive semi-definite.\n";
          std::vector<double> sdiag(old.numRows()+old.numCols());
          std::vector<double> e(old.numRows()+old.numCols());
          std::vector<double> work(old.numRows());
          Matrix<T> u(old.numRows(),old.numCols());
          Matrix<T> s(old.numRows(),old.numCols(),T(0));
          Matrix<T> v(old.numRows(),old.numCols());
          int info = dsvdc ( &old.m_data[0], old.numRows(), old.numRows(), old.numCols(),
                             &sdiag[0], &e[0], &u.m_data[0], old.numRows(), &v.m_data[0],
                             old.numCols(), &work[0], 11 );
          if(info > 0) msg_Error()<<"Error in linpack svd decomposition.\n";
          for(size_t i=0; i<s.numCols(); i++) s(i,i) = abs(sdiag[i])>T(1e-6) ? T(1)/sdiag[i] : 0;
          v *= s*Conjugate(u);
          // msg_Out()<<v<<"\n"<<s<<"\n"<<Transpose(u)<<"\n";
          // exit(1);
          m_data = v.m_data;
          break;
        }
      default: msg_Error()<<"Unexpected Error code in syminv.\n"; break;
      }
      return *this;
    }
    
    /// Calculate the inverse of the symmetric matrix with lower triangular matrix equal to this one.
    Matrix<T> inverseLowerTriangular() {
      return InverseLowerTriangular(*this);
    }  
    // \}
    
    /// \name Element access.
    // \{
    /// Access to element in row \p r, column \p c.
    const T& operator()(size_t r, size_t c) const {return element(r,c);}
    
    /// Access to element in row \p r, column \p c (non-const).
    T& operator()(size_t r, size_t c) {return element(r,c);}
    
    /// Access to element in row \p r, column \p c.
    inline const T& element(size_t r, size_t c) const  {return m_data[r*m_dim + c];}
    
    /// Access to element in row \p r, column \p c (non-const).
    inline T& element(size_t r, size_t c) {return m_data[r*m_dim + c];}
    
    /// Direct access to underlying valarray.
    inline const std::valarray<T>& data() const { return m_data; }
    
    /// Direct access to underlying valarray (non-const).
    inline std::valarray<T>& data() { return m_data; }  
    // \}
    
    
    
    /// \name Misc.
    // \{
    /// Return the number of columns.
    inline size_t numCols() const { return m_dim; }
    
    /// Return the number of rows.
    inline size_t numRows() const { return numCols()==0 ? 0 : numElements()/numCols(); }
    
    /// Return the number of elements in the matrix.
    inline size_t numElements() const { return m_data.size(); } 
    
    /// Determine wheter the matrix is a square matrix.
    inline bool isSquare() const { return numRows()==numCols(); }
    
    /// Determine whether this matrix and \p other have the same number of rows and columns.
    inline bool sameDim(const Matrix<T>& other) const {
      return numRows()==other.numRows() && numCols()==other.numCols();
    }
    
    /// Return a Matrix containing the real part of every element of the matrix.
    Matrix<double> real() const {
      std::valarray<double> data(m_data.size());
      for(size_t i=0; i<data.size(); i++) data[i] = m_data[i].real();
      return Matrix<double>(data, m_dim);
    }

    /// Return a Matrix containing the imaginary part of every elemnet of the matrix.
    Matrix<double> imag() const {
      std::valarray<double> data(m_data.size());
      for(size_t i=0; i<data.size(); i++) data[i] = m_data[i].imag();
      return Matrix<double>(data, m_dim);
    }
    
  
    ///  Sanitize small values to 0.
    Matrix<T> setFuzzyZeroToZero(T eps=T(1e-12)) const {
      Matrix<T> ret = *this;
      return ret.setFuzzyZeroToZeroInline(eps);
    }
    
    /// Sanitize small values to 0 (non-const).
    Matrix<T>& setFuzzyZeroToZeroInline(T eps=T(1e-12)) {
      for(T& v: m_data) {
        if(abs(v)<abs(eps)) v = T(0.);
      }
      return *this;
    }
    
    /// Print the matrix to \p stream.
    void print(std::ostream& stream) const {
      for(size_t r=0; r<numRows(); r++) {
        for(size_t c=0; c<numCols(); c++) {
          stream << element(r,c) << " ";
        }
        stream << "\n";
      }
    }
    // \}
    
  private:
    
    /// \name Internal.
    // \{
    
    /// Does the actual work in corresponding constructors.
    void _init(const std::vector<std::vector<T>>& mat) {
      if(!mat.empty()) {
        m_data.resize(mat.size()*mat.at(0).size());
        m_dim = mat.at(0).size();
        try{
          for(size_t r=0; r<mat.size(); r++) {
            for(size_t c=0; c<m_dim; c++) {
              element(r,c) = mat.at(r).at(c);
            }
          }
        }
        catch(const std::out_of_range& e) {
          THROW(fatal_error, "Tried to initalize matrix from non-compatible structure!");
        }
      }
      else {
        m_data = std::valarray<T>();
        m_dim = 0;
      }
    }
    
    /// Does the actual work in corresponding constructors.
    void _init(const std::vector<std::vector<T>>& mat, size_t rows, size_t cols, T def=T()) {
      m_dim = cols;
      m_data.resize(rows*cols);
      try{
        for(size_t r=0; r<rows; r++) {
          for(size_t c=0; c<cols; c++) {
            element(r,c) = (r<mat.size() && c<mat.at(r).size()) ? mat.at(r).at(c) : def;
          }
        }
      }
      catch(const std::out_of_range& e) {
        THROW(fatal_error, "Tried to initalize matrix from non-compatible structure!");
      }
    }
    
    /// Internal helper, place lower triangular matrix in m_data
    /** Fills the first m_dim*(m_dim+1)/2 entries of m_data with
        element(0,0)
        element(1,0) element(1,1)
        element(2,0) element(2,1) element(2,2)
        ...
        Matrix is not valid after that, this is for use with syminv 
        (and where appropriate). Functions that call this should usually also call
        _fromLowTriang or replace m_data in some other way...
    */
    void _toLowTriang() {
      size_t k = 1;
      for(size_t i=1; i<m_dim; i++) {
        for(size_t j=0; j<=i; j++) {
          std::swap(m_data[k], element(i,j));
          k++;
        }
      }
    }
    
    /// Internal helper, creates m_data from lower triangular matrix \p new_data
    /** Assumes new_data contains values corresponding to
        element(0,0)
        element(1,0) element(1,1)
        element(2,0) element(2,1) element(2,2)
        ...
        and fills m_data correspondingly. Matrix is valid after that.
    */
    void _fromLowTriang(std::valarray<T> new_data) {
      size_t k = 0;
      for(size_t i=0; i<m_dim; i++) {
        for(int j=0; j<=i; j++) {
          std::swap(element(i,j), new_data[k]);
          if(i!=j) element(j,i) =  element(i,j);
          k++;
        }
      }
    }
    // \}
    
  public:
    /// Stream output.
    friend std::ostream& operator<<(std::ostream& stream, const Matrix<T>& matrix) {
      matrix.print(stream);
      return stream;
    }
    
    friend inline Matrix<double> real(const Matrix<T>& matrix) {
      return matrix.real();
    }
    
    friend inline Matrix<double> imag(const Matrix<T>& matrix) {
      return matrix.imag();
    }
    
    /// Calculate the transpose of \p matrix.
    friend inline Matrix<T> Transpose(Matrix<T> matrix) {
      return matrix.transposeInPlace();
    }
    
    
    /// Calculate the trace of \p matrix.
    friend inline T Trace(const Matrix<T>& matrix) {
      return matrix.trace();
    }
    
    /// Calculate the trace of the product of \p a and \p b.
    friend inline T Trace(const Matrix<T>& a, Matrix<T> b) {
      if(a.isSquare() && b.isSquare() && a.sameDim(b)) {
        return (a.m_data*b.transposeInPlace().m_data).sum();
      }
      else{
        return Trace(a*b);
      }
    }

    /// Calculate the trace of the product of \p a and \p b.
    template<template <class X> class ItAble>
    friend inline ItAble<T> Trace(const ItAble<Matrix<T>>& a, Matrix<T> b) {
      ItAble<T> ret(a.size());
      b.transposeInPlace();
      for(size_t i=0; i<a.size(); i++) {
        if(a[i].isSquare() && b.isSquare() && a[i].sameDim(b)) {
          ret[i] = (a[i].m_data*b.m_data).sum();
        }
        else{
          ret[i] = Trace(a[i]*Transpose(b));
        }
      }
      return ret;
    }

    
    
    /// Calculate the conjugate of \p matrix.
    friend inline Matrix<T> Conjugate(Matrix<T> matrix) {
      return matrix.conjugateInPlace();
    }
    
    /// Calculate the hermitian transpose of \p matrix.
    friend inline Matrix<T> HermitianTranspose(Matrix<T> matrix) {
      return matrix.hermitianTransposeInPlace();
    }
    
    /// Calculate matrix exponential of \p x.
    friend Matrix<T> exp(Matrix<T> x) {
      return x.exponentiateInPlace();
    } 
    
    /// Calculate the inverse of the symmetric matrix with lower triangular matrix equal \p matrix.
    friend Matrix<T> InverseLowerTriangular(Matrix<T> matrix) {
      return matrix.inverseLowerTriangularInPlace();
    }
    
    /// Determine whether \p a and \p b have the same number of rows and columns
    friend inline bool sameDim(const Matrix<T>& a, const Matrix<T>& b) {
      return a.sameDim(b);
    }
  };
  
  
  /// Matrix addition.
  template <typename T>
  inline Matrix<T> operator+(Matrix<T> lhs, const Matrix<T>& rhs)
  {
    lhs += rhs;
    return lhs;
  }
  
  
  /// Matrix subtraction.
  template <typename T>
  inline Matrix<T> operator-(Matrix<T> lhs, const Matrix<T>& rhs)
  {
    lhs -= rhs;
    return lhs;
  }
  
  /// Matrix multiplication.
  template <typename T>
  inline Matrix<T> operator*(Matrix<T> lhs, const Matrix<T>& rhs)
  {
    lhs *= rhs;
    return lhs;
  }

  /// Multiplication of matrix with scalar.
  template <typename T, typename Y>
  inline Matrix<T> operator*(Matrix<T> lhs, Y rhs)
  {
    lhs *= rhs;
    return lhs;
  }
  
  /// Multiplication of scalar with matrix.
  template <typename T, typename Y>
  inline Matrix<T> operator*(Y lhs, Matrix<T> rhs)
  {
    rhs *= lhs;
    return rhs;
  }
  
  /// Division of matrix by scalar.
  template <typename T>
  inline Matrix<T> operator/(Matrix<T> lhs, T rhs)
  {
    lhs /= rhs;
    return lhs;
  }
  
  /// Determine whether matrices are equal elementwise.
  template <typename T>
  inline bool operator==(const Matrix<T>& lhs, const Matrix<T>& rhs) {
    return lhs.sameDim(rhs) && (lhs.data()==rhs.data()).min();
  }
  
  /// Determine whether matrices have differing elements.
  template <typename T>
  inline bool operator!=(const Matrix<T>& lhs, const Matrix<T>& rhs) {
    return !(lhs==rhs);
  }
  
  /// Swap algorithm.
  template <typename T>
  void swap(Matrix<T>& a, Matrix<T>& b) {
    a.swap(b);
  }
  
  // usefull typedefs
  /// Typedef for matrix with real elements.
  typedef RESUM::Matrix<double> MatrixD;
  /// Typedef for matrix with complex elements.
  typedef RESUM::Matrix<std::complex<double>> MatrixC;
  
}

/// Swap algorithm.
namespace std {
  /// Specialisation of std::swap.
  template <typename T>
  void swap(RESUM::Matrix<T>& a, RESUM::Matrix<T>& b) {
    a.swap(b);
  }
}

#endif
