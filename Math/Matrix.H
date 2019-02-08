#ifndef MATRIX_H
#define MATRIX_H
/** \file */
#include "ATOOLS/Org/Exception.H"
#include <valarray>
#include <cassert>
#include <complex>


using std::valarray;


/// Basic matrix class.
/**
    Class implementing basic functionality. If you need more, use some
    linear algebra library.
*/
  
template <typename T>
class Matrix {

protected:
  /// Valarray with underlying data.
  valarray<T> m_data;

  /// Number of columns.
  size_t m_dim;
  
public:
  /// \name Standard constructors.
  // \{
  /// Default constructor.
  Matrix() = default;

  /// Copy constructor.
  Matrix<T>(const Matrix<T>& other) = default;

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
  Matrix(size_t r, size_t c, T def=T()) : m_data(r*c), m_dim(c) {}

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
  Matrix(const std::vector<std::vector<T> >& mat, size_t rows, size_t cols, T def=T()) {
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
      m_data = valarray<T>();
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
  Matrix(const std::vector<std::vector<T> >* mat, size_t rows, size_t cols, T def=T()) {
    if(mat) {
      _init(*mat, rows, cols, def);
    }
    else {
      m_data = valarray<T>(def,rows*cols);
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
  static Matrix<T> diagonal(const T& diag, size_t dim, const T& offdiag=T()) {
    Matrix<T> ret(dim,dim,offdiag);
    for(size_t i=0; i<dim;i++) ret(i,i) = diag;
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
  Matrix<T>& operator+=(const Matrix<T>& rhs)
  {
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
    size_t rows = numRows();
    size_t cols = numCols();
    size_t ocols = other.numCols();
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
  Matrix<T>& operator*=(T b)
  {
    m_data *= b;
    return *this;
  }

  /// Division by scalar.
  Matrix<T>& operator/=(T b)
  {
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
  inline const valarray<T>& data() const { return m_data; }
  
  /// Direct access to underlying valarray (non-const).
  inline valarray<T>& data() { return m_data; }  
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
      catch(std::out_of_range) {
        THROW(fatal_error, "Tried to initalize matrix from non-compatible structure!");
        }
    }
    else {
      m_data = valarray<T>();
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
    catch(std::out_of_range) {
      THROW(fatal_error, "Tried to initalize matrix from non-compatible structure!");
    }
  }

public:
  /// Stream output.
  friend std::ostream& operator<<(std::ostream& stream, const Matrix<T>& matrix) {
      matrix.print(stream);
      return stream;
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
template <typename T>
inline Matrix<T> operator*(Matrix<T> lhs, T rhs)
{
  lhs *= rhs;
  return lhs;
}

/// Multiplication of scalar with matrix.
template <typename T>
inline Matrix<T> operator*(T lhs, Matrix<T> rhs)
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

namespace std {
  /// Specialisation of std::swap.
  template <typename T>
  void swap(Matrix<T>& a, Matrix<T>& b) {
    a.swap(b);
  }
}


// usefull typedefs
/// Typedef for matrix with real elements.
typedef Matrix<double> MatrixD;

/// Specialisation of Matrix for complex matrices.
class MatrixC: public Matrix<std::complex<double> > {
public:
  /// \name Standard constructors.
  // \{
  using Matrix::Matrix;
  // \}

  /// \name Special functions for complex matrices.
  // \{
  /// Return a MatrixD containing the real part of every element of the matrix.
  MatrixD real() {
    valarray<double> data(m_data.size());
    for(size_t i=0; i<data.size(); i++) data[i] = m_data[i].real();
    return MatrixD(data, m_dim);
  }

  /// Return a MatrixD containing the imaginary part of every elemnet of the matrix.
  MatrixD imag() {
    valarray<double> data(m_data.size());
    for(size_t i=0; i<data.size(); i++) data[i] = m_data[i].imag();
    return MatrixD(data, m_dim);
  }
  // \}
};

#endif
