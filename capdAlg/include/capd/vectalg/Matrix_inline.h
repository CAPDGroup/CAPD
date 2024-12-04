/// @addtogroup vectalg
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file Matrix_inline.h
///
/// @author The CAPD group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file is part of the CAPD library.  This library is free software;
// you can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation;
// either version 2 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this software; see the file "license.txt".  If not, write to the
// Free Software Foundation, Inc., 59

#ifndef _CAPD_VECTALG_MATRIX_INLINE_H_
#define _CAPD_VECTALG_MATRIX_INLINE_H_

#include "capd/vectalg/Matrix.h"
#include "capd/basicalg/TypeTraits.h"
#include <utility>

namespace capd{
namespace vectalg{


// --------------------- multiplication matrix*vector, matrix*matrix ----------------

template<typename Scalar, __size_type rows,__size_type cols>
inline Vector<Scalar,rows> operator*(const Matrix<Scalar,rows,cols>& m, const Vector<Scalar,cols>& v){
  return matrixByVector< Vector<Scalar,rows> > (m,v);
}

template<typename Scalar, __size_type rows,__size_type cols>
inline Vector<Scalar,rows> operator*(const Matrix<Scalar,rows,cols>& a, const ColumnVector<Scalar,cols>&v){
  return matrixByVector< Vector<Scalar,rows> > (a,v);
}

template<typename Scalar, __size_type rows,__size_type cols>
inline Vector<Scalar,rows> operator*(const Matrix<Scalar,rows,cols>& m, const RowVector<Scalar,cols>& u){
  return matrixByVector< Vector<Scalar,rows> >(m,u);
}

template<typename Scalar, __size_type rows, __size_type cols1, __size_type cols2>
Matrix<Scalar,rows,cols2> operator*(const Matrix<Scalar,rows,cols1>& a1, const Matrix<Scalar,cols1,cols2>& a2){
  return matrixByMatrix< Matrix<Scalar,rows,cols2> > (a1,a2);
}

// ---------------------------- abs(matrix) -----------------------------------------

template<typename Scalar,__size_type rows, __size_type cols>
inline Matrix<Scalar,rows,cols> abs(const Matrix<Scalar,rows,cols>& m){
  return absoluteValue< Matrix<Scalar,rows,cols> > (m);
}

// ---------------------------- unary minus -----------------------------------------

template<typename Scalar,__size_type rows, __size_type cols>
inline Matrix<Scalar,rows,cols> operator-(const Matrix<Scalar,rows,cols>& m){
  return unaryMinus< Matrix<Scalar,rows,cols> > (m);
}

// ---------------------------- matrix + matrix -------------------------------------

template<typename Scalar,__size_type rows, __size_type cols>
inline Matrix<Scalar,rows,cols> operator+(const Matrix<Scalar,rows,cols>& m1, const Matrix<Scalar,rows,cols>& m2){
  return addObjects< Matrix<Scalar,rows,cols> > (m1,m2);
}

// ---------------------------- matrix - matrix -------------------------------------

template<typename Scalar,__size_type rows, __size_type cols>
inline Matrix<Scalar,rows,cols> operator-(const Matrix<Scalar,rows,cols>& m1, const Matrix<Scalar,rows,cols>& m2){
  return subtractObjects< Matrix<Scalar,rows,cols> > (m1,m2);
}

// ---------------------------- matrix *,/ scalar -------------------------------------

template<typename Scalar,__size_type rows, __size_type cols>
inline Matrix<Scalar,rows,cols> operator*(const Matrix<Scalar,rows,cols>& m, const Scalar& s){
  return multiplyObjectScalar< Matrix<Scalar,rows,cols> > (m,s);
}

template<typename Scalar,__size_type rows, __size_type cols>
inline Matrix<Scalar,rows,cols> operator*(const Scalar& s, const Matrix<Scalar,rows,cols>& m){
  return multiplyObjectScalar< Matrix<Scalar,rows,cols> > (m,s);
}

template<typename Scalar,__size_type rows, __size_type cols>
inline Matrix<Scalar,rows,cols> operator/(const Matrix<Scalar,rows,cols>& m, const Scalar& s){
  return divideObjectScalar< Matrix<Scalar,rows,cols> > (m,s);
}

// ---------------------------- matrix + scalar -------------------------------------

template<typename Scalar,__size_type rows, __size_type cols>
inline Matrix<Scalar,rows,cols> operator+(const Matrix<Scalar,rows,cols>& m, const Scalar& s){
  return addObjectScalar< Matrix<Scalar,rows,cols> > (m,s);
}

// ---------------------------- matrix - scalar -------------------------------------

template<typename Scalar,__size_type rows, __size_type cols>
inline Matrix<Scalar,rows,cols> operator-(const Matrix<Scalar,rows,cols>& m, const Scalar& s){
  return subtractObjectScalar< Matrix<Scalar,rows,cols> > (m,s);
}

// ------------------------- inequalities - true if hold true on each coord ---------

template<typename Scalar,__size_type rows, __size_type cols>
inline bool operator<(const Matrix<Scalar,rows,cols>& m1, const Matrix<Scalar,rows,cols>& m2){
  return lessThan(m1,m2);
}

template<typename Scalar,__size_type rows, __size_type cols>
inline bool operator>(const Matrix<Scalar,rows,cols>& m1, const Matrix<Scalar,rows,cols>& m2){
  return greaterThan(m1,m2);
}

template<typename Scalar,__size_type rows, __size_type cols>
inline bool operator<=(const Matrix<Scalar,rows,cols>& m1, const Matrix<Scalar,rows,cols>& m2){
  return lessEqual(m1,m2);
}

template<typename Scalar,__size_type rows, __size_type cols>
inline bool operator>=(const Matrix<Scalar,rows,cols>&m1, const Matrix<Scalar,rows,cols>&m2){
  return greaterEqual(m1,m2);
}

template<typename Scalar,__size_type rows, __size_type cols>
inline bool operator==(const Matrix<Scalar,rows,cols> &a1,const Matrix<Scalar,rows,cols> &a2){
  return equal(a1,a2);
}

template<typename Scalar,__size_type rows, __size_type cols>
inline bool operator!=(const Matrix<Scalar,rows,cols> &a1,const Matrix<Scalar,rows,cols> &a2){
  return notEqual(a1,a2);
}

// ------------------ constructors -------------------------------------------------

template<typename Scalar,__size_type rows, __size_type cols>
inline Matrix<Scalar,rows,cols>::Matrix() : ContainerType()
{}

template<typename Scalar,__size_type rows, __size_type cols>
inline Matrix<Scalar,rows,cols>::Matrix(size_type _rows,size_type _cols)
  : ContainerType(_rows,_cols)
{}

template<typename Scalar,__size_type rows, __size_type cols>
inline Matrix<Scalar,rows,cols>::Matrix(size_type _rows,size_type _cols,bool)
  : ContainerType(_rows,_cols,true)
{}

template<typename Scalar,__size_type rows, __size_type cols>
inline Matrix<Scalar,rows,cols>::Matrix(const Matrix& m)
  : ContainerType(m)
{}

template<typename Scalar,__size_type rows, __size_type cols>
inline Matrix<Scalar,rows,cols>::Matrix(const Dimension& dim)
  : ContainerType(dim)
{}

template<typename Scalar,__size_type rows, __size_type cols>
inline Matrix<Scalar,rows,cols>::Matrix(const Dimension& dim,bool)
  : ContainerType(dim, true)
{}


template<typename Scalar, __size_type rows, __size_type cols>
template<typename S,
        typename std::enable_if<std::is_convertible<S, Scalar>::value && !std::is_same<S, Scalar>::value, int>::type
>
Matrix<Scalar,rows,cols>::Matrix(const Matrix<S,rows,cols>& m)
        : ContainerType(m.numberOfRows(),m.numberOfColumns(),true)
{
    std::copy(m.begin(),m.end(),begin());
}

template<typename Scalar,__size_type rows, __size_type cols>
inline Matrix<Scalar,rows,cols>& Matrix<Scalar,rows,cols>::operator=(const Matrix& a){
   ContainerType::operator=(a);
   return *this;
}


//
//template<typename Scalar, __size_type rows,__size_type cols>
//inline Matrix<Scalar,rows,cols>::Matrix(const Matrix&& a_m)
//  : ContainerType(a_m) {
////  : ContainerType(std::forward<ContainerType>(a_m)) {
//}
//
//template<typename Scalar, __size_type rows,__size_type cols>
//inline Matrix<Scalar,rows,cols>& Matrix<Scalar,rows,cols>::operator=(const Matrix&& a) {
//  //ContainerType::operator= ( std::forward<ContainerType>(a));
//  ContainerType::operator= ( static_cast< const ContainerType &&>(a));
//  return *this;
//}

template<typename Scalar, __size_type rows,__size_type cols>
inline Matrix<Scalar,rows,cols>::Matrix(std::initializer_list< std::initializer_list<ScalarType> > l) 
  : ContainerType(l.size(), l.begin()->size(), false) {
    iterator dest = this->begin();
    for(const auto & v: l){
      if(v.size() != this->numberOfColumns())
        throw std::range_error("Matrix constructor: Rows in initializer list do not have equal sizes.");
      dest = std::copy(v.begin(), v.end(), dest);
    }
 }

// ----------------------------- indexing ---------------------------------------------

template<typename Scalar,__size_type rows, __size_type cols>
inline RowVector<Scalar,cols> Matrix<Scalar,rows,cols>::operator[](size_type i) const {
   return RowVector<Scalar,cols>(data+i*numberOfColumns(),numberOfColumns());
}


template<typename Scalar,__size_type rows, __size_type cols>
inline RowVector<Scalar,cols> Matrix<Scalar,rows,cols>::operator()(size_type i) const{
   return RowVector<Scalar,cols>(data + --i*numberOfColumns(),numberOfColumns());
}

template<typename Scalar,__size_type rows, __size_type cols>
inline RowVector<Scalar,cols> Matrix<Scalar,rows,cols>::row(size_type i) const{
   return RowVector<Scalar,cols>(data+i*numberOfColumns(),numberOfColumns());
}

template<typename Scalar,__size_type rows, __size_type cols>
inline ColumnVector<Scalar,rows> Matrix<Scalar,rows,cols>::column(size_type i) const{
   return ColumnVector<Scalar,rows>(data+i,numberOfColumns(),numberOfRows());
}

template<typename Scalar,__size_type rows, __size_type cols>
inline Scalar& Matrix<Scalar,rows,cols>::operator()(size_type i,size_type j){
   return data[ --i*numberOfColumns() + --j];
}

template<typename Scalar,__size_type rows, __size_type cols>
inline const Scalar& Matrix<Scalar,rows,cols>::operator()(size_type i,size_type j) const{
   return data[ --i*numberOfColumns() + --j];
}

template<typename Scalar,__size_type rows, __size_type cols>
inline Scalar* Matrix<Scalar,rows,cols>::at(size_type i,size_type j){
   return data + --i*numberOfColumns() + --j;
}

template<typename Scalar,__size_type rows, __size_type cols>
inline const Scalar* Matrix<Scalar,rows,cols>::at(size_type i,size_type j) const{
   return data + --i*numberOfColumns() + --j;
}

// ------------------------ iterators ----------------------------- //

template<typename Scalar,__size_type rows, __size_type cols>
inline MatrixIterator< Matrix<Scalar,rows,cols> >
Matrix<Scalar,rows,cols>::beginMatrix() {
  return MatrixIterator< Matrix<Scalar,rows,cols> >(*this,at(1,1));
}

template<typename Scalar,__size_type rows, __size_type cols>
inline MatrixIterator< Matrix<Scalar,rows,cols> >
Matrix<Scalar,rows,cols>::endMatrix(){
  return MatrixIterator< Matrix<Scalar,rows,cols> >(*this,at(numberOfRows(),numberOfColumns()+1));
}

template<typename Scalar,__size_type rows, __size_type cols>
inline MatrixIterator< Matrix<Scalar,rows,cols> > Matrix<Scalar,rows,cols>::beginOfRow(size_type i){
  return MatrixIterator< Matrix<Scalar,rows,cols> >(*this,at(i,1));
}

template<typename Scalar,__size_type rows, __size_type cols>
inline MatrixIterator< Matrix<Scalar,rows,cols> > Matrix<Scalar,rows,cols>::beginOfColumn(size_type j){
  return MatrixIterator< Matrix<Scalar,rows,cols> >(*this,at(1,j));
}

template<typename Scalar,__size_type rows, __size_type cols>
inline MatrixIterator< Matrix<Scalar,rows,cols> > Matrix<Scalar,rows,cols>::endOfRow(size_type i){
  return MatrixIterator< Matrix<Scalar,rows,cols> >(*this,at(i,numberOfColumns()+1));
}

template<typename Scalar,__size_type rows, __size_type cols>
inline MatrixIterator< Matrix<Scalar,rows,cols> > Matrix<Scalar,rows,cols>::endOfColumn(size_type j){
  return MatrixIterator< Matrix<Scalar,rows,cols> >(*this,at(numberOfRows(),j+1));
}


template<typename Scalar,__size_type rows, __size_type cols>
inline const_MatrixIterator< Matrix<Scalar,rows,cols> >
Matrix<Scalar,rows,cols>::beginMatrix() const{
  return const_MatrixIterator< Matrix<Scalar,rows,cols> >(*this,at(1,1));
}

template<typename Scalar,__size_type rows, __size_type cols>
inline const_MatrixIterator< Matrix<Scalar,rows,cols> >
Matrix<Scalar,rows,cols>::endMatrix() const{
  return const_MatrixIterator< Matrix<Scalar,rows,cols> >(*this,at(numberOfRows(),numberOfColumns()+1));
}

template<typename Scalar,__size_type rows, __size_type cols>
inline const_MatrixIterator< Matrix<Scalar,rows,cols> >
Matrix<Scalar,rows,cols>::beginOfRow(size_type i) const{
  return const_MatrixIterator< Matrix<Scalar,rows,cols> >(*this,at(i,1));
}

template<typename Scalar,__size_type rows, __size_type cols>
inline const_MatrixIterator< Matrix<Scalar,rows,cols> >
Matrix<Scalar,rows,cols>::beginOfColumn(size_type j) const{
  return const_MatrixIterator< Matrix<Scalar,rows,cols> >(*this,at(1,j));
}

template<typename Scalar,__size_type rows, __size_type cols>
inline const_MatrixIterator< Matrix<Scalar,rows,cols> >
Matrix<Scalar,rows,cols>::endOfRow(size_type i) const{
  return const_MatrixIterator< Matrix<Scalar,rows,cols> >(*this,at(i,numberOfColumns()+1));
}

template<typename Scalar,__size_type rows, __size_type cols>
inline const_MatrixIterator< Matrix<Scalar,rows,cols> >
Matrix<Scalar,rows,cols>::endOfColumn(size_type j) const{
  return const_MatrixIterator< Matrix<Scalar,rows,cols> >(*this,at(numberOfRows(),j+1));
}

//------------------------ assignments - Scalars ------------------------//

template<typename Scalar,__size_type rows, __size_type cols>
inline Matrix<Scalar,rows,cols>& Matrix<Scalar,rows,cols>::operator=(const Scalar &s){
  return assignFromScalar(*this,s);
}


template<typename Scalar,__size_type rows, __size_type cols>
inline Matrix<Scalar,rows,cols>& Matrix<Scalar,rows,cols>::operator+=(const Scalar &s){
  return addAssignObjectScalar(*this,s);
}

template<typename Scalar,__size_type rows, __size_type cols>
inline Matrix<Scalar,rows,cols>& Matrix<Scalar,rows,cols>::operator-=(const Scalar &s){
  return subtractAssignObjectScalar(*this,s);
}

template<typename Scalar,__size_type rows, __size_type cols>
inline Matrix<Scalar,rows,cols>& Matrix<Scalar,rows,cols>::operator*=(const Scalar &s){
  return multiplyAssignObjectScalar(*this,s);
}

template<typename Scalar,__size_type rows, __size_type cols>
inline Matrix<Scalar,rows,cols>& Matrix<Scalar,rows,cols>::operator/=(const Scalar &s){
  return divideAssignObjectScalar(*this,s);
}

//------------------------- assignments - matrices --------------------//

template<typename Scalar,__size_type rows, __size_type cols>
inline Matrix<Scalar,rows,cols>& Matrix<Scalar,rows,cols>::operator+=(const Matrix& a){
  return addAssignObjectObject(*this,a);
}

template<typename Scalar,__size_type rows, __size_type cols>
inline Matrix<Scalar,rows,cols>& Matrix<Scalar,rows,cols>::operator-=(const Matrix& a){
  return subtractAssignObjectObject(*this,a);
}

///////////////////////////////////////////////////////////////////////////
// trace
///
///  Computes Trace of a given matrix.
///
///  Trace is a sum of all elements on the main diagonal
///
///////////////////////////////////////////////////////////////////////////
template<typename Scalar,__size_type rows, __size_type cols>
typename Matrix<Scalar,cols,rows>::ScalarType trace(const Matrix<Scalar,cols,rows>& A) {

  typedef typename Matrix<Scalar,cols,rows>::ScalarType ScalarType ;
  typedef typename Matrix<Scalar,cols,rows>::size_type size_type ;
  ScalarType T = capd::TypeTraits<ScalarType>::zero();

  for(size_type i = 1; i<= A.numberOfRows(); ++i){
    T += A(i,i);
  }
  return T;
}

///////////////////////////////////////////////////////////////////////////
// secondTrace
///
///  It returns a sum of determinants of all 2x2 matrices
///  from main diagonal
///
///////////////////////////////////////////////////////////////////////////
template<typename Scalar,__size_type rows, __size_type cols>
typename Matrix<Scalar,cols,rows>::ScalarType secondTrace(const Matrix<Scalar,cols,rows> & A) {

  typedef typename Matrix<Scalar,cols,rows>::ScalarType ScalarType ;
  typedef typename Matrix<Scalar,cols,rows>::size_type size_type ;
  ScalarType T = capd::TypeTraits<ScalarType>::zero();

  for(size_type i = 1; i<=A.numberOfRows(); ++i)
    for(size_type j = i+1; j<=A.numberOfColumns(); ++j)
      T += A(i,i)*A(j,j)-A(i,j)*A(j,i);

   return T;
}


}} // namespace capd::vectalg

#endif // _CAPD_VECTALG_MATRIX_INLINE_H_

/// @}
