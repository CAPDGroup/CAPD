/// @addtogroup vectalg
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file MatrixSlice.h
///
/// This class represents a matrix without own container and data
/// it is used for operations on submatrices
///
/// @author Marian Mrozek
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
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


#ifndef _CAPD_VECTALG_MATRIXSLICE_H_ 
#define _CAPD_VECTALG_MATRIXSLICE_H_ 

#include "capd/vectalg/MatrixIterator.h"

/// This class represents a matrix without own container and data
/// it is used for operations on submatrices
template <typename matrix>
class MatrixSlice{
public:
  MatrixSlice(matrix& m,int A_mFirstRow,int A_nLastRow,int A_nFirstCol, int A_nLastCol);
  typedef typename matrix::ScalarType ScalarType;
  typedef typename matrix::size_type size_type;

  ScalarType* at(size_type i,size_type j);
  const ScalarType* at(size_type i,size_type j) const;
  size_type numberOfRows() const;
  size_type numberOfColumns() const;

  size_type rowStride() const;       // difference of pointers when moving to the next row in the same column
  size_type columnStride() const;    // difference of pointers when moving to the next column in the same row

  // iterator selections
  MatrixIterator< MatrixSlice<matrix> > beginMatrix();
  MatrixIterator< MatrixSlice<matrix> > endMatrix();
  MatrixIterator< MatrixSlice<matrix> > beginOfRow(size_type i);
  MatrixIterator< MatrixSlice<matrix> > beginOfColumn(size_type j);
  MatrixIterator< MatrixSlice<matrix> > endOfRow(size_type i);
  MatrixIterator< MatrixSlice<matrix> > endOfColumn(size_type j);

  const_MatrixIterator< MatrixSlice<matrix> > beginMatrix() const;
  const_MatrixIterator< MatrixSlice<matrix> > endMatrix() const;
  const_MatrixIterator< MatrixSlice<matrix> > beginOfRow(size_type i) const;
  const_MatrixIterator< MatrixSlice<matrix> > beginOfColumn(size_type j) const;
  const_MatrixIterator< MatrixSlice<matrix> > endOfRow(size_type i) const;
  const_MatrixIterator< MatrixSlice<matrix> > endOfColumn(size_type j) const;

private:
  matrix& m_baseMatrix;
  size_type m_nFirstRow,m_nLastRow,m_nFirstCol,m_nLastCol;
  //MatrixSlice(MatrixSlice&);
  void operator=(MatrixSlice&);
};


template<typename matrix>
MatrixSlice<matrix>::MatrixSlice(matrix& A_m,int A_nFirstRow,int A_nLastRow,int A_nFirstCol, int A_nLastCol):
  m_baseMatrix(A_m),
  m_nFirstRow(A_nFirstRow),
  m_nLastRow(A_nLastRow),
  m_nFirstCol(A_nFirstCol),
  m_nLastCol(A_nLastCol){
}

template<typename matrix>
inline typename matrix::ScalarType* MatrixSlice<matrix>::at(size_type i,size_type j)
{
   return m_baseMatrix.at(m_nFirstRow+i-1,m_nFirstCol+j-1);
}

template<typename matrix>
inline const typename matrix::ScalarType* MatrixSlice<matrix>::at(size_type i,size_type j) const
{
   return m_baseMatrix.at(m_nFirstRow+i-1,m_nFirstCol+j-1);
}

template<typename matrix>
inline typename MatrixSlice<matrix>::size_type MatrixSlice<matrix>::numberOfRows() const{
  return m_nLastRow-m_nFirstRow+1;
}

template<typename matrix>
inline typename MatrixSlice<matrix>::size_type MatrixSlice<matrix>::numberOfColumns() const{
  return m_nLastCol-m_nFirstCol+1;
}

template<typename matrix>
inline typename MatrixSlice<matrix>::size_type MatrixSlice<matrix>::rowStride() const{
  return m_baseMatrix.rowStride();
}

template<typename matrix>
inline typename MatrixSlice<matrix>::size_type MatrixSlice<matrix>::columnStride() const{
  return m_baseMatrix.columnStride();
}


template<typename matrix>
inline MatrixIterator< MatrixSlice<matrix> > MatrixSlice<matrix>::beginMatrix(){
  return MatrixIterator< MatrixSlice<matrix> >(*this,at(1,1));
}
template<typename matrix>
inline MatrixIterator< MatrixSlice<matrix> > MatrixSlice<matrix>::endMatrix(){
  return MatrixIterator< MatrixSlice<matrix> >(*this,at(numberOfRows(),numberOfColumns()+1));
}
template<typename matrix>
inline MatrixIterator< MatrixSlice<matrix> > MatrixSlice<matrix>::beginOfRow(size_type i){
  return MatrixIterator< MatrixSlice<matrix> >(*this,at(i,1));
}
template<typename matrix>
inline MatrixIterator< MatrixSlice<matrix> > MatrixSlice<matrix>::beginOfColumn(size_type j){
  return MatrixIterator< MatrixSlice<matrix> >(*this,at(1,j));
}
template<typename matrix>
inline MatrixIterator< MatrixSlice<matrix> > MatrixSlice<matrix>::endOfRow(size_type i){
  return MatrixIterator< MatrixSlice<matrix> >(*this,at(i,numberOfColumns()+1));
}
template<typename matrix>
inline MatrixIterator< MatrixSlice<matrix> > MatrixSlice<matrix>::endOfColumn(size_type j){
  return MatrixIterator< MatrixSlice<matrix> >(*this,at(numberOfRows(),j+1));
}

template<typename matrix>
inline const_MatrixIterator< MatrixSlice<matrix> > MatrixSlice<matrix>::beginMatrix() const{
  return const_MatrixIterator< MatrixSlice<matrix> >(*this,at(1,1));
}
template<typename matrix>
inline const_MatrixIterator< MatrixSlice<matrix> > MatrixSlice<matrix>::endMatrix() const{
  return const_MatrixIterator< MatrixSlice<matrix> >(*this,at(numberOfRows(),numberOfColumns()+1));
}
template<typename matrix>
inline const_MatrixIterator< MatrixSlice<matrix> > MatrixSlice<matrix>::beginOfRow(size_type i) const{
  return const_MatrixIterator< MatrixSlice<matrix> >(*this,at(i,1));
}
template<typename matrix>
inline const_MatrixIterator< MatrixSlice<matrix> > MatrixSlice<matrix>::beginOfColumn(size_type j) const{
  return const_MatrixIterator< MatrixSlice<matrix> >(*this,at(1,j));
}
template<typename matrix>
inline const_MatrixIterator< MatrixSlice<matrix> > MatrixSlice<matrix>::endOfRow(size_type i) const{
  return const_MatrixIterator< MatrixSlice<matrix> >(*this,at(i,numberOfColumns()+1));
}
template<typename matrix>
inline const_MatrixIterator< MatrixSlice<matrix> > MatrixSlice<matrix>::endOfColumn(size_type j) const{
  return const_MatrixIterator< MatrixSlice<matrix> >(*this,at(numberOfRows(),j+1));
}

#endif // _CAPD_VECTALG_MATRIXSLICE_H_ 

/// @}
