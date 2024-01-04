/// @addtogroup vectalg
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file Matrix.hpp
///
/// @author Marian Mrozek, Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_VECTALG_MATRIX_HPP_
#define _CAPD_VECTALG_MATRIX_HPP_

#include <vector>
#include <stack>
#include <sstream>
#include <algorithm>

#include "capd/basicalg/minmax.h"
#include "capd/vectalg/Matrix.h"

#include "capd/vectalg/RowVector.hpp"
#include "capd/vectalg/ColumnVector.hpp"
#include "capd/vectalg/Vector.hpp"
#include "capd/basicalg/TypeTraits.h"

#include "capd/matrixAlgorithms/floatMatrixAlgorithms.hpp"
#include "capd/matrixAlgorithms/intMatrixAlgorithms.hpp"
#include "capd/vectalg/algebraicOperations.hpp"

namespace capd{
namespace vectalg{


//-------------------------- constructors ------------------------//

template<typename Scalar, __size_type rows, __size_type cols>
Matrix<Scalar,rows,cols>::Matrix(const MatrixSlice<Matrix> &m)
   : ContainerType(m.numberOfRows(),m.numberOfColumns(),true)
{
  if(numberOfRows() && numberOfColumns()){
   for (size_type i=1; i<=numberOfRows(); ++i)
    for (size_type j=1; j<=numberOfColumns(); ++j)
      (*this)(i,j) = *m.at(i,j);
  }
}

template<typename Scalar, __size_type rows, __size_type cols>
Matrix<Scalar,rows,cols>::Matrix(size_type _rows,size_type _cols,const Scalar dane[]) : ContainerType(_rows,_cols,true)
{
  std::copy(dane,dane+this->size(),begin());
}

template<typename Scalar, __size_type rows, __size_type cols>
Matrix<Scalar,rows,cols>::Matrix(const Scalar dane[]) : ContainerType() {
  std::copy(dane,dane+this->size(),begin());
}

template<typename Scalar, __size_type rows, __size_type cols>
Matrix<Scalar,rows,cols>::Matrix(const char data[]) : ContainerType()
{
   std::istringstream str(data);
   str >> *this;
}

template<typename Scalar, __size_type rows, __size_type cols>
Matrix<Scalar,rows,cols>::Matrix(const std::string & data) : ContainerType()
{
   std::istringstream str(data);
   str >> *this;
}

template<typename Scalar, __size_type rows, __size_type cols>
template<__size_type dataRows,__size_type dataCols>
Matrix<Scalar,rows,cols>::Matrix(const Scalar (&data)[dataRows][dataCols]): ContainerType(dataRows, dataCols, true)
{
  for (__size_type i = 0; i < dataRows; ++i) {
    for (__size_type j = 0; j < dataCols; ++j) {
      (*this)[i][j] = data[i][j];
    }
  }
}

template<typename Scalar,__size_type rows, __size_type cols>
inline void setDimensionInArray(Matrix<Scalar,rows,cols>* , __size_type , __size_type , __size_type ){
}

template<typename Scalar>
inline void setDimensionInArray(Matrix<Scalar,0,0>* t, __size_type N, __size_type r, __size_type c){
  for(__size_type i=0;i<N;++i) t[i].resize(r,c);
}

template<typename Scalar,__size_type rows, __size_type cols>
Matrix<Scalar,rows,cols>* Matrix<Scalar,rows,cols>::makeArray(size_type N, size_type r, size_type c)
{
  Matrix* result = new Matrix[N];
  setDimensionInArray(result,N,r,c);
  return result;
}

//------------- operations on matrices -------------------------//

template<typename Scalar, __size_type rows, __size_type cols>
Matrix<Scalar,rows,cols> Matrix<Scalar,rows,cols>::Identity(size_type dim)
{
   // the matrix is filled by zeros
   Matrix<Scalar,rows,cols> temp(dim,dim);
   for(size_type i=1;i<=dim;++i)
      temp(i,i) = TypeTraits<ScalarType>::one();
   return temp;
}


template<typename Scalar, __size_type rows, __size_type cols>
void Matrix<Scalar,rows,cols>::setToIdentity()
{
   if(numberOfRows()!=numberOfColumns())
      throw std::range_error("Matrix<Scalar,rows,cols>::setToIdentity: rows!=cols");
   this->clear();
   for(size_type i=1;i<=numberOfRows();++i)
      (*this)(i,i) = TypeTraits<ScalarType>::one();
}


template<typename Scalar, __size_type rows, __size_type cols>
void Matrix<Scalar,rows,cols>::transpose()
{
   if(numberOfRows()!=numberOfColumns())
      throw std::runtime_error("Cannot call x.transpose() for nonsquare matrix");

   for(size_type i=1;i<=numberOfRows();++i)
   for(size_type j=i+1;j<=numberOfRows();++j)
   {
      Scalar a = (*this)(i,j);
      (*this)(i,j) = (*this)(j,i);
      (*this)(j,i) = a;
   }
}


template<typename Scalar, __size_type rows, __size_type cols>
Matrix<Scalar,cols,rows> transpose(const Matrix<Scalar,rows,cols> &a)
{
  typedef typename Matrix<Scalar,cols,rows>::size_type size_type;
  Matrix<Scalar,cols,rows> temp(a.numberOfColumns(),a.numberOfRows(),true);
   for(size_type i=1;i<=a.numberOfColumns();++i)
      for(size_type j=1;j<=a.numberOfRows();++j)
         temp(i,j)= a(j,i);
   return temp;
}

template<typename Scalar, __size_type rows, __size_type cols>
inline Matrix<Scalar,cols,rows> Transpose(const Matrix<Scalar,rows,cols> &a){
  return transpose(a);
}
// ----------------- input - output ------------------------------ //

template<typename Scalar, __size_type rows, __size_type cols>
std::ostream &operator<<(std::ostream &out,const Matrix<Scalar,rows,cols> &a)
{
  typedef typename Matrix<Scalar,cols,rows>::size_type size_type;
   out << "{";
   if(a.numberOfColumns()>0){
     if(a.numberOfRows()>1) out << std::endl;
     if(a.numberOfRows()>0) out << a[0];
     for(size_type i=1;i<a.numberOfRows();++i)
     {
        out << "," << std::endl << a[i];
     }
     if(a.numberOfRows()>1) out << std::endl;
   }
   out << "}";
   return out;
}


template<typename Scalar, __size_type rows, __size_type cols>
std::istream &operator>>(std::istream &inp, Matrix<Scalar,rows,cols> &a){
  typedef typename Matrix<Scalar,rows,cols>::RowVectorType VectorType;

  std::stack<VectorType> st;
  VectorType v(cols);
  int ch=0;
  typename Matrix<Scalar,rows,cols>::size_type dimension=0;
  while('{'!=(ch=inp.get()) && ch!=EOF);
  if(ch!=EOF){
    while(' '==(ch=inp.peek()) || ch=='\n') inp.get();
    if(ch!='}'){ // otherwise we have an ampty matrix
      inp >> v;
      st.push(v);
      dimension=v.dimension();
      do{
        do{
          ch=inp.get();
        }while(ch==' ' || ch=='\t' || ch=='\n');
        if(ch==','){
          inp >> v;
          if(dimension!=v.dimension())
            throw std::ios_base::failure("Rows of incompatible dimensions found when reading a matrix");
          st.push(v);
        }
      }while(ch!='}' && ch!=EOF);
    }
  }

  int n = static_cast<int>(st.size());
  a.resize(n,dimension);
  for(int i=n-1;i>=0;--i)
  {
    v=st.top();
    a[i]=v;
    st.pop();
  }

  if(inp.eof())
    throw std::ios_base::failure("EOF encountered when reading a matrix");
  return inp;
}


template<typename Scalar, __size_type rows,__size_type cols>
std::string cppReprezentation(const Matrix<Scalar,cols,rows> & A, const std::string& varName,
			      const std::string& typeName)
{
  std::stringstream out;
  out << "capd::vectalg::Matrix<" << typeName << ", " << rows << ", " << cols << "> " << varName << "(";

  if (A.numberOfRows() > 0 && A.numberOfColumns() > 0) {
    out << "(" << typeName << "[" << A.numberOfRows() << "][" << A.numberOfColumns() << "])" << A;
  } else {
    out << A.numberOfRows() << ", " << A.numberOfColumns();
  }

  out << ");";
  std::string str = out.str();
  std::replace(str.begin(), str.end(), '\n', ' ');

  return str;
}

}} // namespace capd::vectalg

#endif // _CAPD_VECTALG_MATRIX_HPP_

/// @}
