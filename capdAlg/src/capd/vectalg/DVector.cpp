
/////////////////////////////////////////////////////////////////////////////
/// @file vectalg/DVector.cpp
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include <iostream>
#include "capd/vectalg/Dimension.h"
#include "capd/vectalg/Vector.hpp"
#include "capd/vectalg/RowVector.hpp"
#include "capd/vectalg/ColumnVector.hpp"

namespace capd{
  namespace vectalg{

template class Vector<double,CAPD_DEFAULT_DIMENSION>;

template Vector<double,CAPD_DEFAULT_DIMENSION> abs<double,CAPD_DEFAULT_DIMENSION> (const Vector<double,CAPD_DEFAULT_DIMENSION> &v);
template Vector<double,CAPD_DEFAULT_DIMENSION> operator- <double,CAPD_DEFAULT_DIMENSION> (const Vector<double,CAPD_DEFAULT_DIMENSION> &v);
template Vector<double,CAPD_DEFAULT_DIMENSION> operator+ <double,CAPD_DEFAULT_DIMENSION> (const Vector<double,CAPD_DEFAULT_DIMENSION> &v1,const Vector<double,CAPD_DEFAULT_DIMENSION> &v2);
template Vector<double,CAPD_DEFAULT_DIMENSION> operator- <double,CAPD_DEFAULT_DIMENSION> (const Vector<double,CAPD_DEFAULT_DIMENSION> &v1,const Vector<double,CAPD_DEFAULT_DIMENSION> &v2);
template double operator* <double,CAPD_DEFAULT_DIMENSION> (const Vector<double,CAPD_DEFAULT_DIMENSION> &v1,const Vector<double,CAPD_DEFAULT_DIMENSION> &v2);

template Vector<double,CAPD_DEFAULT_DIMENSION> operator* <double,double,CAPD_DEFAULT_DIMENSION> (const Vector<double,CAPD_DEFAULT_DIMENSION> &v,const double &s);
template Vector<double,CAPD_DEFAULT_DIMENSION> operator* <double,double,CAPD_DEFAULT_DIMENSION> (const double &s,const Vector<double,CAPD_DEFAULT_DIMENSION> &v);
template Vector<double,CAPD_DEFAULT_DIMENSION> operator/ <double,double,CAPD_DEFAULT_DIMENSION> (const Vector<double,CAPD_DEFAULT_DIMENSION> &v, const double &s);

template Vector<double,CAPD_DEFAULT_DIMENSION> operator* <double,long,CAPD_DEFAULT_DIMENSION>(const Vector<double,CAPD_DEFAULT_DIMENSION> &v,const long &s);
template Vector<double,CAPD_DEFAULT_DIMENSION> operator* <double,long,CAPD_DEFAULT_DIMENSION>(const long &s,const Vector<double,CAPD_DEFAULT_DIMENSION> &v);
template Vector<double,CAPD_DEFAULT_DIMENSION> operator/ <double,long,CAPD_DEFAULT_DIMENSION>(const Vector<double,CAPD_DEFAULT_DIMENSION> &v, const long &s);

template Vector<double,CAPD_DEFAULT_DIMENSION> operator* <double,int,CAPD_DEFAULT_DIMENSION>(const Vector<double,CAPD_DEFAULT_DIMENSION> &v,const int &s);
template Vector<double,CAPD_DEFAULT_DIMENSION> operator* <double,int,CAPD_DEFAULT_DIMENSION>(const int &s,const Vector<double,CAPD_DEFAULT_DIMENSION> &v);
template Vector<double,CAPD_DEFAULT_DIMENSION> operator/ <double,int,CAPD_DEFAULT_DIMENSION>(const Vector<double,CAPD_DEFAULT_DIMENSION> &v, const int &s);

template Vector<double,CAPD_DEFAULT_DIMENSION> operator+ <double,CAPD_DEFAULT_DIMENSION> (const Vector<double,CAPD_DEFAULT_DIMENSION> &v,const double &s);
template Vector<double,CAPD_DEFAULT_DIMENSION> operator- <double,CAPD_DEFAULT_DIMENSION> (const Vector<double,CAPD_DEFAULT_DIMENSION> &v,const double &s);
template bool operator < <double,CAPD_DEFAULT_DIMENSION> (const Vector<double,CAPD_DEFAULT_DIMENSION> &v1,const Vector<double,CAPD_DEFAULT_DIMENSION> &v2);
template bool operator > <double,CAPD_DEFAULT_DIMENSION> (const Vector<double,CAPD_DEFAULT_DIMENSION> &v1,const Vector<double,CAPD_DEFAULT_DIMENSION> &v2);
template bool operator<= <double,CAPD_DEFAULT_DIMENSION> (const Vector<double,CAPD_DEFAULT_DIMENSION> &v1,const Vector<double,CAPD_DEFAULT_DIMENSION> &v2);
template bool operator>= <double,CAPD_DEFAULT_DIMENSION> (const Vector<double,CAPD_DEFAULT_DIMENSION> &v1,const Vector<double,CAPD_DEFAULT_DIMENSION> &v2);
template std::ostream &operator<< <double,CAPD_DEFAULT_DIMENSION> (std::ostream &out, const Vector<double,CAPD_DEFAULT_DIMENSION> &v);
template std::istream &operator>> <double,CAPD_DEFAULT_DIMENSION>(std::istream &inp, Vector<double,CAPD_DEFAULT_DIMENSION> &v);

template void subtractObjects<>(const Vector<double,CAPD_DEFAULT_DIMENSION>& v1,const Vector<double,CAPD_DEFAULT_DIMENSION>& v2, Vector<double,CAPD_DEFAULT_DIMENSION>& result);
template void addObjects<>(const Vector<double,CAPD_DEFAULT_DIMENSION>& v1,const Vector<double,CAPD_DEFAULT_DIMENSION>& v2, Vector<double,CAPD_DEFAULT_DIMENSION>& result);

template RowVector<double, CAPD_DEFAULT_DIMENSION>& assignObjectObject<>(RowVector<double, CAPD_DEFAULT_DIMENSION>&, RowVector<double, CAPD_DEFAULT_DIMENSION> const&);
template RowVector<long double, CAPD_DEFAULT_DIMENSION>& assignObjectObject<>(RowVector<long double, CAPD_DEFAULT_DIMENSION>&, RowVector<long double, CAPD_DEFAULT_DIMENSION> const&);
template ColumnVector<double, CAPD_DEFAULT_DIMENSION>& assignObjectObject<>(ColumnVector<double, CAPD_DEFAULT_DIMENSION>&, ColumnVector<double, CAPD_DEFAULT_DIMENSION> const&);
template ColumnVector<double, CAPD_DEFAULT_DIMENSION>& assignObjectObject<>(ColumnVector<double, CAPD_DEFAULT_DIMENSION>&, Vector<double, CAPD_DEFAULT_DIMENSION> const&);
template ColumnVector<long double, CAPD_DEFAULT_DIMENSION>& assignObjectObject<>(ColumnVector<long double, CAPD_DEFAULT_DIMENSION>&, ColumnVector<long double, CAPD_DEFAULT_DIMENSION> const&);
template ColumnVector<long double, CAPD_DEFAULT_DIMENSION>& assignObjectObject<>(ColumnVector<long double, CAPD_DEFAULT_DIMENSION>&, Vector<long double, CAPD_DEFAULT_DIMENSION> const&);
template Vector<double, CAPD_DEFAULT_DIMENSION>& assignObjectObject<>(Vector<double, CAPD_DEFAULT_DIMENSION>&, Vector<double, CAPD_DEFAULT_DIMENSION> const&);

template Vector<double, CAPD_DEFAULT_DIMENSION> multiplyObjectScalar<>(RowVector<double, CAPD_DEFAULT_DIMENSION> const&, double const&);
template Vector<long double, CAPD_DEFAULT_DIMENSION> multiplyObjectScalar<>(RowVector<long double, CAPD_DEFAULT_DIMENSION> const&, long double const&);
template RowVector<double, CAPD_DEFAULT_DIMENSION>& divideAssignObjectScalar<>(RowVector<double, CAPD_DEFAULT_DIMENSION>&, double const&);
template RowVector<long double, CAPD_DEFAULT_DIMENSION>& divideAssignObjectScalar<>(RowVector<long double, CAPD_DEFAULT_DIMENSION>&, long double const&);
template Vector<double, CAPD_DEFAULT_DIMENSION> divideObjectScalar<Vector<double, CAPD_DEFAULT_DIMENSION>, Vector<double, CAPD_DEFAULT_DIMENSION>, long>(Vector<double, CAPD_DEFAULT_DIMENSION> const&, long const&);
template Vector<long double, CAPD_DEFAULT_DIMENSION> divideObjectScalar<Vector<long double, CAPD_DEFAULT_DIMENSION>, Vector<long double, CAPD_DEFAULT_DIMENSION>, long>(Vector<long double, CAPD_DEFAULT_DIMENSION> const&, long const&);

template Vector<double, CAPD_DEFAULT_DIMENSION> unaryMinus<>(ColumnVector<double, CAPD_DEFAULT_DIMENSION> const&);
template Vector<long double, CAPD_DEFAULT_DIMENSION> unaryMinus<>(ColumnVector<long double, CAPD_DEFAULT_DIMENSION> const&);
template std::ostream & printVector<>(std::ostream & str, const ColumnVector<double, CAPD_DEFAULT_DIMENSION> &, int, int);
template std::ostream & printVector<>(std::ostream & str, const ColumnVector<long double, CAPD_DEFAULT_DIMENSION> &, int, int);

template std::string cppReprezentation<>(Vector<double, CAPD_DEFAULT_DIMENSION> const&, std::string const&, std::string const&);
template std::string cppReprezentation<>(Vector<long double, CAPD_DEFAULT_DIMENSION> const&, std::string const&, std::string const&);

template Vector<double, CAPD_DEFAULT_DIMENSION> multiplyObjectScalar<>(ColumnVector<double, CAPD_DEFAULT_DIMENSION> const&, double const&);
template Vector<double, CAPD_DEFAULT_DIMENSION> multiplyObjectScalar<>(Vector<double, CAPD_DEFAULT_DIMENSION>  const&, double const&);
template Vector<double, CAPD_DEFAULT_DIMENSION> multiplyObjectScalar<>(Vector<double, CAPD_DEFAULT_DIMENSION>  const&, long const&);
template Vector<long double, CAPD_DEFAULT_DIMENSION> multiplyObjectScalar<>(Vector<long double, CAPD_DEFAULT_DIMENSION>  const&, long double const&);
template Vector<long double, CAPD_DEFAULT_DIMENSION> multiplyObjectScalar<>(Vector<long double, CAPD_DEFAULT_DIMENSION>  const&, long  const&);
template Vector<long double, CAPD_DEFAULT_DIMENSION> multiplyObjectScalar<>(ColumnVector<long double, CAPD_DEFAULT_DIMENSION>  const&, long double const&);
template ColumnVector<double, CAPD_DEFAULT_DIMENSION>& divideAssignObjectScalar<>(ColumnVector<double, CAPD_DEFAULT_DIMENSION>&, double const&);
template ColumnVector<long double, CAPD_DEFAULT_DIMENSION>& divideAssignObjectScalar<>(ColumnVector<long double, CAPD_DEFAULT_DIMENSION>&, long double const&);
template void clear<>(ColumnVector<double, CAPD_DEFAULT_DIMENSION>&);
template void clear<>(ColumnVector<long double, CAPD_DEFAULT_DIMENSION>&);
template Vector<double, CAPD_DEFAULT_DIMENSION>& multiplyAssignObjectScalarAddObject<>(Vector<double, CAPD_DEFAULT_DIMENSION>&, double const&, Vector<double, CAPD_DEFAULT_DIMENSION> const&);
template Vector<long double, CAPD_DEFAULT_DIMENSION>& multiplyAssignObjectScalarAddObject<>(Vector<long double, CAPD_DEFAULT_DIMENSION>&, long double const&, Vector<long double, CAPD_DEFAULT_DIMENSION> const&);
template RowVector<double, CAPD_DEFAULT_DIMENSION>& subtractAssignObjectObject<>(RowVector<double, CAPD_DEFAULT_DIMENSION>&, Vector<double, CAPD_DEFAULT_DIMENSION> const&);
template Vector<double, CAPD_DEFAULT_DIMENSION>& addAssignObjectObject<>(Vector<double, CAPD_DEFAULT_DIMENSION>&, Vector<double, CAPD_DEFAULT_DIMENSION> const&);

typedef Vector<double,CAPD_DEFAULT_DIMENSION> DVector;
template DVector divideObjectScalar<>(const DVector & v1, const double & s);
}}  // end of namespace capd::vectalg
