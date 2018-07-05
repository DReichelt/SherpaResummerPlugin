#ifndef MATRIX_H
#define MATRIX_H

#include <valarray>
#include <cassert>
#include <iostream>

using std::valarray;
using std::slice;

template <typename T>
class Matrix {
protected:
  std::valarray<T> m_data;
  size_t m_dim;
public:
  // constructors  
  // need that to keep stl containers happy...
  Matrix() {}

  // directly provide underlying data and number of cols
  Matrix(const std::valarray<T> d, size_t c) : m_data(d), m_dim(c) {
    assert(d.size()%c == 0);
  }
  // empty matrix with r rows and c cols
  Matrix(size_t r, size_t c) : m_data(r*c), m_dim(c) {}
  // build matrix from nested vectors, determine number of rows and cols from vector
  Matrix(const std::vector<std::vector<T> >& mat) : m_dim(mat.size()), m_data(mat.size()*mat.size()) {
    try{
      for(size_t i=0; i<m_dim; i++) {
        for(size_t j=0; j<m_dim; j++) {
          element(i,j) = mat.at(i).at(j);
        }
      }
    }
    catch(std::out_of_range) {
      THROW(fatal_error, "Tried to initalize matrix from non-compatible structure!")
    }
  }
  // build matrix from nested vectors, but fill with zeros to match dim rows and cols
  Matrix(const std::vector<std::vector<T> >& mat, size_t dim, T def) : m_dim(dim), m_data(dim*dim) {
    try{
      for(size_t i=0; i<m_dim; i++) {
        for(size_t j=0; j<m_dim; j++) {
          element(i,j) = (i<mat.size() && j<mat.size()) ? mat.at(i).at(j) : def;
        }
      }
    }
    catch(std::out_of_range) {
      THROW(fatal_error, "Tried to initalize matrix from non-compatible structure!")
    }
  }

  static Matrix<T> diagonal(const T& diag, const T& offdiag, size_t dim) {
    Matrix<T> ret(dim,dim);
    for(size_t i=0;i<dim;++i) {
      for(size_t j=0;j<dim;j++) {
        ret(i,j) = i==j?diag:offdiag;
      }
    }
    return ret;
  }
  
  // same as above but from pointer to nested vector
  Matrix(const std::vector<std::vector<T> >* mat) {
    if(mat) {
      m_dim = mat->size();
      m_data = valarray<T>(m_dim*m_dim);
      try{
        for(size_t i=0; i<m_dim; i++) {
          for(size_t j=0; j<m_dim; j++) {
            element(i,j) = mat->at(i).at(j);
          }
        }
      }
      catch(std::out_of_range) {
        THROW(fatal_error, "Tried to initalize matrix from non-compatible structure!")
          }
    }
    else {
      m_dim = 0;
      m_data = valarray<T>();
    }
  }
  // same as above but from pointer to nested vector
  Matrix(const std::vector<std::vector<T> >* mat, size_t dim, T def) : m_dim(dim), m_data(dim*dim){
    if(mat) {
      try{
        for(size_t i=0; i<m_dim; i++) {
          for(size_t j=0; j<m_dim; j++) {
            element(i,j) = (i<mat->size() && j<mat->size()) ? mat->at(i).at(j) : def;
          }
        }
      }
      catch(std::out_of_range) {
        THROW(fatal_error, "Tried to initalize matrix from non-compatible structure!")
          }
    }
    else {
      m_dim = dim;
      m_data = valarray<T>(def, dim*dim);
    }
  }
  

  // basic operations
  Matrix<T>& operator+=(const Matrix<T>& rhs)
  {
    assert(m_dim == rhs.dim());
    m_data += rhs.m_data;
    return *this;
  }

  Matrix<T>& operator-=(const Matrix<T>& rhs) {
    assert(m_dim == rhs.dim());
    m_data -= rhs.m_data;
    return *this;
  }

  Matrix<T> operator*=(const Matrix<T>& other) {
    assert(m_dim == other.dim());
    valarray<T> tmp = m_data;
    for(size_t i = 0; i < m_dim; ++i)
      for(size_t j = 0; j < m_dim; ++j) {
        const valarray<T>& row = tmp[slice(m_dim*i, m_dim, 1)];
        const valarray<T>& col = other.m_data[slice(j, m_dim, m_dim)];
        element(i,j) = (row * col).sum();
      }
    return *this;
  }

  
  Matrix<T>& operator*=(T b)
  {
    m_data *= b;
    return *this;
  }

  Matrix<T>& operator/=(T b)
  {
    m_data /= b;
    return *this;
  }

  
  // accessors
  T& operator()(size_t r, size_t c) {return element(r,c);}
  const T& operator()(size_t r, size_t c) const {return element(r,c);}
  inline T& element(size_t r, size_t c) {return m_data[r*m_dim + c];}
  inline const T& element(size_t r, size_t c) const  {return m_data[r*m_dim + c];}


  
  // properties of the matrix
  T trace() const {
    return m_data[std::slice(0, m_dim, m_dim+1)].sum();
  }


  // misc
  inline size_t dim() const { return m_dim; }
  inline size_t numElements() const { return m_dim*m_dim; } 

  void print(std::ostream& stream) const {
    for(size_t i=0; i<m_dim; i++) {
      for(size_t j=0; j<m_dim; j++) {
        stream << element(i,j) << " ";
      }
      stream << "\n";
    }
  }

  
  
  friend std::ostream& operator<< (std::ostream& stream, const Matrix<T>& matrix) {
      matrix.print(stream);
      return stream;
  }
  

};

template <typename T>
inline Matrix<T> operator+(Matrix<T> lhs, const Matrix<T>& rhs)
{
  lhs += rhs;
  return lhs;
}

template <typename T>
inline Matrix<T> operator-(Matrix<T> lhs, const Matrix<T>& rhs)
{
  lhs -= rhs;
  return lhs;
}

template <typename T>
inline Matrix<T> operator*(Matrix<T> lhs, const Matrix<T>& rhs)
{
  lhs *= rhs;
  return lhs;
}

template <typename T>
inline Matrix<T> operator*(Matrix<T> lhs, T rhs)
{
  lhs *= rhs;
  return lhs;
}

template <typename T>
inline Matrix<T> operator*(T lhs, Matrix<T> rhs)
{
  rhs *= lhs;
  return rhs;
}


template <typename T>
inline Matrix<T> operator/(Matrix<T> lhs, T rhs)
{
  lhs /= rhs;
  return lhs;
}



// usefull typedef
typedef Matrix<double> MatrixD;

// specialisation for complex matrices
class MatrixC: public Matrix<std::complex<double> > {
public:
  using Matrix::Matrix;
  MatrixD real() {
    valarray<double> data(m_data.size());
    for(int i=0; i<data.size(); i++) data[i] = m_data[i].real();
    return MatrixD(data, m_dim);
  }
  MatrixD imag() {
    valarray<double> data(m_data.size());
    for(int i=0; i<data.size(); i++) data[i] = m_data[i].imag();
    return MatrixD(data, m_dim);
  }
};



#endif
