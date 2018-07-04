#include "Tools/Hard_Matrix.H"
#include "ATOOLS/Org/Message.H"

using namespace PHASIC;

namespace RESUM {
  std::ostream &operator<<(std::ostream &s,const Hard_Matrix &h)
  {
    s<<"Hard_Matrix(this="<<&h<<"){\n";
    for (size_t i(0);i<h.m_id.size();++i)
    s<<"  c_"<<i<<" = "<<h.m_id[i]<<"\n";
    for (size_t i(0);i<h.size();++i)
    for (size_t j(0);j<h[i].size();++j)
    s<<"  h["<<i<<"]["<<j<<"] = "<<h[i][j]<<"\n";
    return s<<"}";
  }

  int size( const Hard_Matrix &h)
  {
    int result = h.size();
    return result;
  }

  int bas( const Hard_Matrix &h , const int i, const int k)
  {
    int result = h.m_id[i][k];
    return result;
  }

  double access( const Hard_Matrix &h , const int i , const int j )
  {
  return double(real(h[i][j]));
  }

}
