/// @addtogroup vectalg
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file Matrix.h
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

#ifndef _CAPD_VECTALG_MATRIX_H_
#define _CAPD_VECTALG_MATRIX_H_

#include <iostream>
#include "capd/vectalg/Vector.h"
#include "capd/vectalg/RowVector.h"
#include "capd/vectalg/ColumnVector.h"
#include "capd/vectalg/MatrixContainer.h"
#include "capd/vectalg/MatrixIterator.h"
#include "capd/vectalg/MatrixSlice.h"
#include "capd/settings/compilerSetting.h"


namespace capd{
namespace vectalg{

template<typename Scalar, __size_type rows, __size_type cols>
class Matrix;

template<typename Scalar, __size_type rows,__size_type cols1, __size_type cols2>
Matrix<Scalar,rows,cols2> operator*(const Matrix<Scalar,rows,cols1>&, const Matrix<Scalar,cols1,cols2>&);

template<typename Scalar, __size_type rows,__size_type cols>
Matrix<Scalar,cols,rows> transpose(const Matrix<Scalar,rows,cols>&);

//Deprecated
template<typename Scalar, __size_type rows,__size_type cols>
Matrix<Scalar,cols,rows> Transpose(const Matrix<Scalar,rows,cols>&);

template<typename Scalar, __size_type rows, __size_type cols>
std::ostream &operator<<(std::ostream&, const Matrix<Scalar,rows,cols>&);

template<typename Scalar, __size_type rows,__size_type cols>
std::istream &operator>>(std::istream&, Matrix<Scalar,rows,cols>&);

// ########################################################################### //

template<typename Scalar, __size_type rows, __size_type cols>
class Matrix : public MatrixContainer<Scalar,rows,cols>
{
public:
  typedef Scalar ScalarType;
  typedef MatrixContainer<Scalar,rows,cols> ContainerType;
  typedef typename ContainerType::iterator iterator;
  typedef typename ContainerType::const_iterator const_iterator;
  typedef typename ContainerType::Dimension Dimension;
  typedef __size_type size_type;
  typedef __difference_type difference_type;

  typedef Vector<Scalar,cols> RowVectorType;
  typedef Vector<Scalar,rows> ColumnVectorType;
  typedef RowVector<Scalar,cols> RefRowVectorType;
  typedef ColumnVector<Scalar,rows> RefColumnVectorType;
  typedef Matrix<Scalar,rows,cols> MatrixType;

  //constructors
  Matrix();
  Matrix(size_type _rows,size_type _cols);               //assigns 0 to every coordinate
  Matrix(const ScalarType data[]);
  Matrix(size_type _rows,size_type _cols, const ScalarType data[]);
  explicit Matrix(const char data[]);
  explicit Matrix(const std::string & data);
  Matrix(const Matrix& m);                     // copying constructor
  Matrix(const MatrixSlice<MatrixType>& m);    // copying from matrix slice


  template<typename S,
          typename std::enable_if<std::is_convertible<S, ScalarType>::value && !std::is_same<S, ScalarType>::value, int>::type = 0 >
  Matrix(const Matrix<S,rows,cols> &);

  Matrix(const Dimension& dim);           //assigns 0 to every coordinate
  Matrix(const Dimension& d,bool);
  Matrix(size_type _rows,size_type _cols,bool);

  template<__size_type dataRows,__size_type dataCols>
  Matrix(const Scalar (&data)[dataRows][dataCols]);
  Matrix& operator= (Matrix&& a)  noexcept = default;         //move a matrix
  Matrix(Matrix&& m) noexcept = default;                     // move constructor
  Matrix(std::initializer_list< std::initializer_list<ScalarType> > l);


  //assignments - matrices
  Matrix& operator= (const Matrix& a);         //assign a matrix
  Matrix& operator+=(const Matrix& a);         //increase by a matrix
  Matrix& operator-=(const Matrix& a);         //decrease by a matrix

  //assignments - Scalars
  Matrix& operator= (const Scalar& s);         //assign a Scalar
  Matrix& operator+=(const Scalar& s);         //increase by a Scalar
  Matrix& operator-=(const Scalar& s);         //decrease by a Scalar
  Matrix& operator*=(const Scalar& s);         //rescale by multiplying
  Matrix& operator/=(const Scalar& s);         //rescale by dividing

  // iterator selections
  MatrixIterator< MatrixType > beginMatrix();
  MatrixIterator< MatrixType > endMatrix();
  MatrixIterator< MatrixType > beginOfRow(size_type i);
  MatrixIterator< MatrixType > beginOfColumn(size_type j);
  MatrixIterator< MatrixType > endOfRow(size_type i);
  MatrixIterator< MatrixType > endOfColumn(size_type j);

  const_MatrixIterator< MatrixType > beginMatrix() const;
  const_MatrixIterator< MatrixType > endMatrix() const;
  const_MatrixIterator< MatrixType > beginOfRow(size_type i) const;
  const_MatrixIterator< MatrixType > beginOfColumn(size_type j) const;
  const_MatrixIterator< MatrixType > endOfRow(size_type i) const;
  const_MatrixIterator< MatrixType > endOfColumn(size_type j) const;

  using ContainerType::begin;
  using ContainerType::end;

  //indexing
  RowVector<Scalar,cols> operator[](size_type i) const; // i-th row as a vector
  RowVector<Scalar,cols> operator()(size_type i) const; // (i-1)-th row as a vector
  Scalar& operator()(size_type i,size_type j);                 //returns reference to the [i-1][j-1] entry
  const Scalar& operator()(size_type i,size_type j) const;     //returns const reference to the [i-1][j-1] entry
  Scalar* at(size_type i,size_type j);                         //returns pointer to the [i-1][j-1] entry
  const Scalar* at(size_type i,size_type j) const;             //returns pointer to the [i-1][j-1] entry

  //operations on matrices
  static Matrix Identity(size_type dim);              // returns the identity matrix
  void setToIdentity();                               // if square matrix, changes it to the identity matrix
  RowVector<Scalar,cols> row(size_type i) const;       // i-th row as a vector
  ColumnVector<Scalar,rows> column(size_type j) const; // returns j-th column
  void transpose();

  difference_type rowStride()const {
    return numberOfColumns();
  } // difference of pointers when moving to the next row in the same column
  difference_type columnStride()const{
    return 1;
  }  // difference of pointers when moving to the next column in the same row

  template<typename U>
  struct rebind {
      typedef Matrix<U,rows,cols> other;
  };

  using ContainerType::numberOfRows;
  using ContainerType::numberOfColumns;
  using ContainerType::dimension;

  friend std::istream &operator>> <> (std::istream &inp, MatrixType &a);
  static Matrix* makeArray(size_type N, size_type r, size_type c);

// Deprecated
   void Transpose() { transpose(); }
   using ContainerType::resize;
   static size_type degree() {return 1;} // required interface for DynSys
protected:

  using ContainerType::data;
}; //the end of class Matrix


///  Computes Trace of a given matrix.
template<typename Scalar, __size_type rows,__size_type cols>
typename Matrix<Scalar,cols,rows>::ScalarType trace(const Matrix<Scalar,cols,rows>& A);

///  It returns a sum of determinants of all 2x2 matrix
template<typename Scalar, __size_type rows,__size_type cols>
typename Matrix<Scalar,cols,rows>::ScalarType secondTrace(const Matrix<Scalar,cols,rows> & A);

/// It serializes a matrix - gives text reprezentation which can be compiled
template<typename Scalar, __size_type rows,__size_type cols>
std::string cppReprezentation(const Matrix<Scalar,cols,rows> & A, const std::string& varName,
			      const std::string& typeName);

}} // namespace capd::vectalg

#include "capd/vectalg/Matrix_inline.h"

#endif // _CAPD_VECTALG_MATRIX_H_

/// @}
