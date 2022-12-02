
/////////////////////////////////////////////////////////////////////////////
/// @file vectalg/DMatrix.cpp
///
/// @author The CAPD Group
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

template class Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>;
template class RowVector<double,CAPD_DEFAULT_DIMENSION>;
template class ColumnVector<double,CAPD_DEFAULT_DIMENSION>;

template Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> abs <double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&);
template Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> operator- <double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&);
template Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> operator+ <double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&);
template Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> operator- <double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&);
template Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> operator* <double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&);
template Vector<double,CAPD_DEFAULT_DIMENSION> operator* <double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const Vector<double,CAPD_DEFAULT_DIMENSION>&);
template Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> operator* <double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const double&);
template Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> operator* <double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const double&, const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&);
template Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> operator/ <double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const double&);
template Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> operator+ <double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const double&);
template Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> operator- <double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const double&);
template bool operator< <double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&);
template bool operator> <double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&);
template bool operator<= <double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&);
template bool operator>= <double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&);
template Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> Transpose <double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&);
template std::ostream &operator<< <double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (std::ostream&, const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&);
template std::istream &operator>> <double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> (std::istream&, Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&);


template void matrixByVector<> (const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const Vector<double,CAPD_DEFAULT_DIMENSION>&,Vector<double,CAPD_DEFAULT_DIMENSION>&);
template void matrixByVector<> (const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const ColumnVector<double,CAPD_DEFAULT_DIMENSION>&,Vector<double,CAPD_DEFAULT_DIMENSION>&);
template void matrixByVector<> (const Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const ColumnVector<long double,CAPD_DEFAULT_DIMENSION>&,Vector<long double,CAPD_DEFAULT_DIMENSION>&);
template void matrixByMatrix<> (const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&,Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&);

template void matrixByMatrix<> (const Matrix<int,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const Matrix<int,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&,Matrix<int,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&);
template void matrixByVector<> (const Matrix<int,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const Vector<int,CAPD_DEFAULT_DIMENSION>&,Vector<int,CAPD_DEFAULT_DIMENSION>&);
template void matrixByVector<> (const Matrix<int,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const ColumnVector<int,CAPD_DEFAULT_DIMENSION>&,Vector<int,CAPD_DEFAULT_DIMENSION>&);

template void subtractObjects<>(const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>& v1,const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>& v2, Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>& result);
template void addObjects<>(const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>& v1,const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>& v2, Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>& result);
template Vector<double,CAPD_DEFAULT_DIMENSION>& addAssignMatrixByVector<>(Vector<double,CAPD_DEFAULT_DIMENSION>& u,const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>& A, const Vector<double,CAPD_DEFAULT_DIMENSION>& v);
template Vector<double,CAPD_DEFAULT_DIMENSION>& subtractAssignMatrixByVector<>(Vector<double,CAPD_DEFAULT_DIMENSION>& u,const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>& A, const Vector<double,CAPD_DEFAULT_DIMENSION>& v);

// RowVector

template Vector<double,CAPD_DEFAULT_DIMENSION> operator+<double,CAPD_DEFAULT_DIMENSION>(
      const Vector<double,CAPD_DEFAULT_DIMENSION>&,
      const RowVector<double,CAPD_DEFAULT_DIMENSION>&
   );
template Vector<double,CAPD_DEFAULT_DIMENSION> operator+<double,CAPD_DEFAULT_DIMENSION>(
      const RowVector<double,CAPD_DEFAULT_DIMENSION>&,
      const Vector<double,CAPD_DEFAULT_DIMENSION>&
   );
template Vector<double,CAPD_DEFAULT_DIMENSION> operator+<double,CAPD_DEFAULT_DIMENSION>(
      const RowVector<double,CAPD_DEFAULT_DIMENSION>&,
      const RowVector<double,CAPD_DEFAULT_DIMENSION>&
   );

template Vector<double,CAPD_DEFAULT_DIMENSION> operator-<double,CAPD_DEFAULT_DIMENSION>(
      const Vector<double,CAPD_DEFAULT_DIMENSION>&,
      const RowVector<double,CAPD_DEFAULT_DIMENSION>&
   );
template Vector<double,CAPD_DEFAULT_DIMENSION> operator-<double,CAPD_DEFAULT_DIMENSION>(
      const RowVector<double,CAPD_DEFAULT_DIMENSION>&,
      const Vector<double,CAPD_DEFAULT_DIMENSION>&
   );
template Vector<double,CAPD_DEFAULT_DIMENSION> operator-<double,CAPD_DEFAULT_DIMENSION>(
      const RowVector<double,CAPD_DEFAULT_DIMENSION>&,
      const RowVector<double,CAPD_DEFAULT_DIMENSION>&
   );

template double operator*<double,CAPD_DEFAULT_DIMENSION>(
      const Vector<double,CAPD_DEFAULT_DIMENSION>&,
      const RowVector<double,CAPD_DEFAULT_DIMENSION>&
   );
template double operator*<double,CAPD_DEFAULT_DIMENSION>(
      const RowVector<double,CAPD_DEFAULT_DIMENSION>&,
      const Vector<double,CAPD_DEFAULT_DIMENSION>&
   );
template double operator*<double,CAPD_DEFAULT_DIMENSION>(
      const RowVector<double,CAPD_DEFAULT_DIMENSION>&,
      const RowVector<double,CAPD_DEFAULT_DIMENSION>&
   );

template Vector<double,CAPD_DEFAULT_DIMENSION> operator-<double,CAPD_DEFAULT_DIMENSION>(const RowVector<double,CAPD_DEFAULT_DIMENSION>&);
template Vector<double,CAPD_DEFAULT_DIMENSION> operator*<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>(
      const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, const RowVector<double,CAPD_DEFAULT_DIMENSION>&
   );

template Vector<double,CAPD_DEFAULT_DIMENSION> operator*<double,CAPD_DEFAULT_DIMENSION>(const double&, const RowVector<double,CAPD_DEFAULT_DIMENSION>&);
template Vector<double,CAPD_DEFAULT_DIMENSION> operator*<double,CAPD_DEFAULT_DIMENSION>(const RowVector<double,CAPD_DEFAULT_DIMENSION>&, const double&);
template Vector<double,CAPD_DEFAULT_DIMENSION> operator/<double,CAPD_DEFAULT_DIMENSION>(const RowVector<double,CAPD_DEFAULT_DIMENSION>&, const double&);

template std::ostream& operator<< <double,CAPD_DEFAULT_DIMENSION>(std::ostream&, const RowVector<double,CAPD_DEFAULT_DIMENSION>&);


// ColumnVector

template Vector<double,CAPD_DEFAULT_DIMENSION> operator+<double,CAPD_DEFAULT_DIMENSION>(
      const Vector<double,CAPD_DEFAULT_DIMENSION>&,
      const ColumnVector<double,CAPD_DEFAULT_DIMENSION>&
   );
template Vector<double,CAPD_DEFAULT_DIMENSION> operator+<double,CAPD_DEFAULT_DIMENSION>(
      const ColumnVector<double,CAPD_DEFAULT_DIMENSION>&,
      const Vector<double,CAPD_DEFAULT_DIMENSION>&
   );
template Vector<double,CAPD_DEFAULT_DIMENSION> operator+<double,CAPD_DEFAULT_DIMENSION>(
      const ColumnVector<double,CAPD_DEFAULT_DIMENSION>&,
      const ColumnVector<double,CAPD_DEFAULT_DIMENSION>&
   );

template Vector<double,CAPD_DEFAULT_DIMENSION> operator-<double,CAPD_DEFAULT_DIMENSION>(
      const Vector<double,CAPD_DEFAULT_DIMENSION>&,
      const ColumnVector<double,CAPD_DEFAULT_DIMENSION>&
   );
template Vector<double,CAPD_DEFAULT_DIMENSION> operator-<double,CAPD_DEFAULT_DIMENSION>(
      const ColumnVector<double,CAPD_DEFAULT_DIMENSION>&,
      const Vector<double,CAPD_DEFAULT_DIMENSION>&
   );
template Vector<double,CAPD_DEFAULT_DIMENSION> operator-<double,CAPD_DEFAULT_DIMENSION>(
      const ColumnVector<double,CAPD_DEFAULT_DIMENSION>&,
      const ColumnVector<double,CAPD_DEFAULT_DIMENSION>&
   );

template double operator*<double,CAPD_DEFAULT_DIMENSION>(
      const Vector<double,CAPD_DEFAULT_DIMENSION>&,
      const ColumnVector<double,CAPD_DEFAULT_DIMENSION>&
   );
template double operator*<double,CAPD_DEFAULT_DIMENSION>(
      const ColumnVector<double,CAPD_DEFAULT_DIMENSION>&,
      const Vector<double,CAPD_DEFAULT_DIMENSION>&
   );
template double operator*<double,CAPD_DEFAULT_DIMENSION>(
      const ColumnVector<double,CAPD_DEFAULT_DIMENSION>&,
      const ColumnVector<double,CAPD_DEFAULT_DIMENSION>&
   );

template Vector<double,CAPD_DEFAULT_DIMENSION> operator-<double,CAPD_DEFAULT_DIMENSION>(const ColumnVector<double,CAPD_DEFAULT_DIMENSION>&);
template Vector<double,CAPD_DEFAULT_DIMENSION> operator*<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>(
      const Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&,
      const ColumnVector<double,CAPD_DEFAULT_DIMENSION>&
   );

template Vector<double,CAPD_DEFAULT_DIMENSION> operator*<double,CAPD_DEFAULT_DIMENSION>(const double&, const ColumnVector<double,CAPD_DEFAULT_DIMENSION>&);
template Vector<double,CAPD_DEFAULT_DIMENSION> operator*<double,CAPD_DEFAULT_DIMENSION>(const ColumnVector<double,CAPD_DEFAULT_DIMENSION>&, const double&);
template Vector<double,CAPD_DEFAULT_DIMENSION> operator/<double,CAPD_DEFAULT_DIMENSION>(const ColumnVector<double,CAPD_DEFAULT_DIMENSION>&, const double&);

template std::ostream& operator<< <double,CAPD_DEFAULT_DIMENSION>(std::ostream&, const ColumnVector<double,CAPD_DEFAULT_DIMENSION>&);

template Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>& multiplyAssignObjectScalarAddObject<>(Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, double const&, Matrix<double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> const&);
template Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>& multiplyAssignObjectScalarAddObject<>(Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>&, long double const&, Matrix<long double,CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> const&);

template void Matrix<int, CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>::setToIdentity();
template Matrix<int, CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> Matrix<int, CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>::Identity(size_type);
template Matrix<int, CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION>::Matrix(MatrixSlice<Matrix<int, CAPD_DEFAULT_DIMENSION,CAPD_DEFAULT_DIMENSION> > const&);
template std::string cppReprezentation<>(Matrix<int, CAPD_DEFAULT_DIMENSION, CAPD_DEFAULT_DIMENSION> const&, std::string const&, std::string const&);
template std::string cppReprezentation<>(Matrix<double, CAPD_DEFAULT_DIMENSION, CAPD_DEFAULT_DIMENSION> const&, std::string const&, std::string const&);

  }}
