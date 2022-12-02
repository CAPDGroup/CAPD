
/////////////////////////////////////////////////////////////////////////////
/// @file vectalg/LDMatrix.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include "capd/vectalg/Dimension.h"
#include "capd/vectalg/Vector.hpp"
#include "capd/vectalg/ColumnVector.hpp"
#include "capd/vectalg/RowVector.hpp"
#include "capd/vectalg/Matrix.hpp"

namespace capd{
  namespace vectalg{

template class Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>;
template class RowVector<long double,CAPD_DEFAULT_DIMENSION>;
template class ColumnVector<long double,CAPD_DEFAULT_DIMENSION>;

template Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> abs <long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&);
template Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> operator- <long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&);
template Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> operator+ <long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&);
template Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> operator- <long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&);
template Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> operator* <long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&);
template Vector<long double,CAPD_DEFAULT_DIMENSION> operator* <long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const Vector<long double,CAPD_DEFAULT_DIMENSION>&);
template Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> operator* <long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const long double&);
template Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> operator* <long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const long double&, const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&);
template Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> operator/ <long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const long double&);
template Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> operator+ <long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const long double&);
template Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> operator- <long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const long double&);
template bool operator< <long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&);
template bool operator> <long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&);
template bool operator<= <long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&);
template bool operator>= <long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&);
template Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> Transpose <long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&);
template std::ostream &operator<< <long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (std::ostream&, const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&);
template std::istream &operator>> <long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (std::istream&, Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&);


template void matrixByVector<> (const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const Vector<long double,CAPD_DEFAULT_DIMENSION>&,Vector<long double,CAPD_DEFAULT_DIMENSION>&);
template void matrixByMatrix<> (const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&,Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&);
template void subtractObjects<>(const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>& v1,const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>& v2, Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>& result);
template void addObjects<>(const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>& v1,const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>& v2, Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>& result);
template Vector<long double,CAPD_DEFAULT_DIMENSION>& addAssignMatrixByVector<>(Vector<long double,CAPD_DEFAULT_DIMENSION>& u,const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>& A, const Vector<long double,CAPD_DEFAULT_DIMENSION>& v);
template Vector<long double,CAPD_DEFAULT_DIMENSION>& subtractAssignMatrixByVector<>(Vector<long double,CAPD_DEFAULT_DIMENSION>& u,const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>& A, const Vector<long double,CAPD_DEFAULT_DIMENSION>& v);

// RowVector

template Vector<long double,CAPD_DEFAULT_DIMENSION> operator+<long double,CAPD_DEFAULT_DIMENSION>(
      const Vector<long double,CAPD_DEFAULT_DIMENSION>&,
      const RowVector<long double,CAPD_DEFAULT_DIMENSION>&
   );
template Vector<long double,CAPD_DEFAULT_DIMENSION> operator+<long double,CAPD_DEFAULT_DIMENSION>(
      const RowVector<long double,CAPD_DEFAULT_DIMENSION>&,
      const Vector<long double,CAPD_DEFAULT_DIMENSION>&
   );
template Vector<long double,CAPD_DEFAULT_DIMENSION> operator+<long double,CAPD_DEFAULT_DIMENSION>(
      const RowVector<long double,CAPD_DEFAULT_DIMENSION>&,
      const RowVector<long double,CAPD_DEFAULT_DIMENSION>&
   );

template Vector<long double,CAPD_DEFAULT_DIMENSION> operator-<long double,CAPD_DEFAULT_DIMENSION>(
      const Vector<long double,CAPD_DEFAULT_DIMENSION>&,
      const RowVector<long double,CAPD_DEFAULT_DIMENSION>&
   );
template Vector<long double,CAPD_DEFAULT_DIMENSION> operator-<long double,CAPD_DEFAULT_DIMENSION>(
      const RowVector<long double,CAPD_DEFAULT_DIMENSION>&,
      const Vector<long double,CAPD_DEFAULT_DIMENSION>&
   );
template Vector<long double,CAPD_DEFAULT_DIMENSION> operator-<long double,CAPD_DEFAULT_DIMENSION>(
      const RowVector<long double,CAPD_DEFAULT_DIMENSION>&,
      const RowVector<long double,CAPD_DEFAULT_DIMENSION>&
   );

template long double operator*<long double,CAPD_DEFAULT_DIMENSION>(
      const Vector<long double,CAPD_DEFAULT_DIMENSION>&,
      const RowVector<long double,CAPD_DEFAULT_DIMENSION>&
   );
template long double operator*<long double,CAPD_DEFAULT_DIMENSION>(
      const RowVector<long double,CAPD_DEFAULT_DIMENSION>&,
      const Vector<long double,CAPD_DEFAULT_DIMENSION>&
   );
template long double operator*<long double,CAPD_DEFAULT_DIMENSION>(
      const RowVector<long double,CAPD_DEFAULT_DIMENSION>&,
      const RowVector<long double,CAPD_DEFAULT_DIMENSION>&
   );

template Vector<long double,CAPD_DEFAULT_DIMENSION> operator-<long double,CAPD_DEFAULT_DIMENSION>(const RowVector<long double,CAPD_DEFAULT_DIMENSION>&);
template Vector<long double,CAPD_DEFAULT_DIMENSION> operator*<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>(
      const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const RowVector<long double,CAPD_DEFAULT_DIMENSION>&
   );

template Vector<long double,CAPD_DEFAULT_DIMENSION> operator*<long double,CAPD_DEFAULT_DIMENSION>(const long double&, const RowVector<long double,CAPD_DEFAULT_DIMENSION>&);
template Vector<long double,CAPD_DEFAULT_DIMENSION> operator*<long double,CAPD_DEFAULT_DIMENSION>(const RowVector<long double,CAPD_DEFAULT_DIMENSION>&, const long double&);
template Vector<long double,CAPD_DEFAULT_DIMENSION> operator/<long double,CAPD_DEFAULT_DIMENSION>(const RowVector<long double,CAPD_DEFAULT_DIMENSION>&, const long double&);

template std::ostream& operator<< <long double,CAPD_DEFAULT_DIMENSION>(std::ostream&, const RowVector<long double,CAPD_DEFAULT_DIMENSION>&);


// ColumnVector

template Vector<long double,CAPD_DEFAULT_DIMENSION> operator+<long double,CAPD_DEFAULT_DIMENSION>(
      const Vector<long double,CAPD_DEFAULT_DIMENSION>&,
      const ColumnVector<long double,CAPD_DEFAULT_DIMENSION>&
   );
template Vector<long double,CAPD_DEFAULT_DIMENSION> operator+<long double,CAPD_DEFAULT_DIMENSION>(
      const ColumnVector<long double,CAPD_DEFAULT_DIMENSION>&,
      const Vector<long double,CAPD_DEFAULT_DIMENSION>&
   );
template Vector<long double,CAPD_DEFAULT_DIMENSION> operator+<long double,CAPD_DEFAULT_DIMENSION>(
      const ColumnVector<long double,CAPD_DEFAULT_DIMENSION>&,
      const ColumnVector<long double,CAPD_DEFAULT_DIMENSION>&
   );

template Vector<long double,CAPD_DEFAULT_DIMENSION> operator-<long double,CAPD_DEFAULT_DIMENSION>(
      const Vector<long double,CAPD_DEFAULT_DIMENSION>&,
      const ColumnVector<long double,CAPD_DEFAULT_DIMENSION>&
   );
template Vector<long double,CAPD_DEFAULT_DIMENSION> operator-<long double,CAPD_DEFAULT_DIMENSION>(
      const ColumnVector<long double,CAPD_DEFAULT_DIMENSION>&,
      const Vector<long double,CAPD_DEFAULT_DIMENSION>&
   );
template Vector<long double,CAPD_DEFAULT_DIMENSION> operator-<long double,CAPD_DEFAULT_DIMENSION>(
      const ColumnVector<long double,CAPD_DEFAULT_DIMENSION>&,
      const ColumnVector<long double,CAPD_DEFAULT_DIMENSION>&
   );

template long double operator*<long double,CAPD_DEFAULT_DIMENSION>(
      const Vector<long double,CAPD_DEFAULT_DIMENSION>&,
      const ColumnVector<long double,CAPD_DEFAULT_DIMENSION>&
   );
template long double operator*<long double,CAPD_DEFAULT_DIMENSION>(
      const ColumnVector<long double,CAPD_DEFAULT_DIMENSION>&,
      const Vector<long double,CAPD_DEFAULT_DIMENSION>&
   );
template long double operator*<long double,CAPD_DEFAULT_DIMENSION>(
      const ColumnVector<long double,CAPD_DEFAULT_DIMENSION>&,
      const ColumnVector<long double,CAPD_DEFAULT_DIMENSION>&
   );

template Vector<long double,CAPD_DEFAULT_DIMENSION> operator-<long double,CAPD_DEFAULT_DIMENSION>(const ColumnVector<long double,CAPD_DEFAULT_DIMENSION>&);
template Vector<long double,CAPD_DEFAULT_DIMENSION> operator*<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>(
      const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&,
      const ColumnVector<long double,CAPD_DEFAULT_DIMENSION>&
   );

template Vector<long double,CAPD_DEFAULT_DIMENSION> operator*<long double,CAPD_DEFAULT_DIMENSION>(const long double&, const ColumnVector<long double,CAPD_DEFAULT_DIMENSION>&);
template Vector<long double,CAPD_DEFAULT_DIMENSION> operator*<long double,CAPD_DEFAULT_DIMENSION>(const ColumnVector<long double,CAPD_DEFAULT_DIMENSION>&, const long double&);
template Vector<long double,CAPD_DEFAULT_DIMENSION> operator/<long double,CAPD_DEFAULT_DIMENSION>(const ColumnVector<long double,CAPD_DEFAULT_DIMENSION>&, const long double&);

template std::ostream& operator<< <long double,CAPD_DEFAULT_DIMENSION>(std::ostream&, const ColumnVector<long double,CAPD_DEFAULT_DIMENSION>&);

  }}

