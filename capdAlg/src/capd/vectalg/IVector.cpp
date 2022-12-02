
/////////////////////////////////////////////////////////////////////////////
/// @file vectalg/IVector.cpp
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include <iostream>

#include "capd/intervals/minmax_interval.h"
#include "capd/vectalg/Dimension.h"
#include "capd/vectalg/Vector.hpp"
#include "capd/vectalg/ColumnVector.hpp"
#include "capd/vectalg/RowVector.hpp"
#include "capd/intervals/lib.h"
#include "capd/vectalg/Vector_Interval.hpp"

namespace capd{
  namespace vectalg{


template class Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION>;

template Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> abs<capd::DInterval,CAPD_DEFAULT_DIMENSION> (const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v);
template Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> operator- <capd::DInterval,CAPD_DEFAULT_DIMENSION> (const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v);
template Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> operator+ <capd::DInterval,CAPD_DEFAULT_DIMENSION>(const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v1,const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v2);
template Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> operator- <capd::DInterval,CAPD_DEFAULT_DIMENSION>(const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v1,const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v2);
template capd::DInterval operator* <capd::DInterval,CAPD_DEFAULT_DIMENSION> (const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v1,const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v2);

template Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> operator* <capd::DInterval,capd::DInterval,CAPD_DEFAULT_DIMENSION>(const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v,const capd::DInterval &s);
template Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> operator* <capd::DInterval,capd::DInterval,CAPD_DEFAULT_DIMENSION>(const capd::DInterval &s,const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v);
template Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> operator/ <capd::DInterval,capd::DInterval,CAPD_DEFAULT_DIMENSION>(const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v, const capd::DInterval &s);

template Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> operator* <capd::DInterval,double,CAPD_DEFAULT_DIMENSION>(const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v,const double &s);
template Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> operator* <capd::DInterval,double,CAPD_DEFAULT_DIMENSION>(const double &s,const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v);
template Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> operator/ <capd::DInterval,double,CAPD_DEFAULT_DIMENSION>(const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v, const double &s);

template Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> operator* <capd::DInterval,long,CAPD_DEFAULT_DIMENSION>(const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v,const long &s);
template Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> operator* <capd::DInterval,long,CAPD_DEFAULT_DIMENSION>(const long &s,const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v);
template Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> operator/ <capd::DInterval,long,CAPD_DEFAULT_DIMENSION>(const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v, const long &s);

template Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> operator* <capd::DInterval,int,CAPD_DEFAULT_DIMENSION>(const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v,const int &s);
template Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> operator* <capd::DInterval,int,CAPD_DEFAULT_DIMENSION>(const int &s,const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v);
template Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> operator/ <capd::DInterval,int,CAPD_DEFAULT_DIMENSION>(const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v, const int &s);

template Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> operator+ <capd::DInterval,CAPD_DEFAULT_DIMENSION> (const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v, const capd::DInterval &s);
template Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> operator- <capd::DInterval,CAPD_DEFAULT_DIMENSION>(const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v,const capd::DInterval &s);
template bool operator < <capd::DInterval,CAPD_DEFAULT_DIMENSION> (const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v1,const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v2);
template bool operator > <capd::DInterval,CAPD_DEFAULT_DIMENSION> (const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v1,const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v2);
template bool operator<= <capd::DInterval,CAPD_DEFAULT_DIMENSION> (const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v1,const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v2);
template bool operator>= <capd::DInterval,CAPD_DEFAULT_DIMENSION> (const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v1,const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v2);
template std::ostream &operator<< <capd::DInterval,CAPD_DEFAULT_DIMENSION> (std::ostream &out, const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v);
template std::istream &operator>> <capd::DInterval,CAPD_DEFAULT_DIMENSION>(std::istream &inp, Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> &v);

template void subtractObjects<>(const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION>&,const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION>&, Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION>&);
template void addObjects<>(const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION>&,const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION>&, Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION>&);
template void addObjects<>(const Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION>&,const ColumnVector<capd::DInterval,CAPD_DEFAULT_DIMENSION>&, Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION>&);
template Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> unaryMinus<>(ColumnVector<capd::DInterval,CAPD_DEFAULT_DIMENSION> const&);
template Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION>& addAssignObjectScalar<>(Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION>& v1, const capd::DInterval&);
template Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> addObjectScalar<>(Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> const&, capd::DInterval const&);
template void subtractObjects<>(Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> const&, ColumnVector<capd::DInterval,CAPD_DEFAULT_DIMENSION> const&, Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION>&);
template ColumnVector<capd::DInterval, CAPD_DEFAULT_DIMENSION>& assignObjectObject<>(ColumnVector<capd::DInterval, CAPD_DEFAULT_DIMENSION>&, Vector<capd::DInterval, CAPD_DEFAULT_DIMENSION> const&);
template ColumnVector<capd::DInterval, CAPD_DEFAULT_DIMENSION>& assignObjectObject<>(ColumnVector<capd::DInterval, CAPD_DEFAULT_DIMENSION>&, ColumnVector<capd::DInterval, CAPD_DEFAULT_DIMENSION> const&);
template Vector<capd::DInterval, CAPD_DEFAULT_DIMENSION>& assignObjectObject<>(Vector<capd::DInterval, CAPD_DEFAULT_DIMENSION>&, Vector<double, CAPD_DEFAULT_DIMENSION> const&);
template bool normalize<>(ColumnVector<capd::DInterval,  CAPD_DEFAULT_DIMENSION>&);

template RowVector<capd::DInterval, CAPD_DEFAULT_DIMENSION>& assignObjectObject<>(RowVector<capd::DInterval, CAPD_DEFAULT_DIMENSION>&, Vector<capd::DInterval, CAPD_DEFAULT_DIMENSION> const&);
template RowVector<capd::DInterval, CAPD_DEFAULT_DIMENSION>& assignObjectObject<>(RowVector<capd::DInterval, CAPD_DEFAULT_DIMENSION>&, RowVector<capd::DInterval, CAPD_DEFAULT_DIMENSION> const&);

typedef Vector<capd::DInterval,CAPD_DEFAULT_DIMENSION> IVector;
typedef Vector<capd::DInterval::BoundType,CAPD_DEFAULT_DIMENSION> DVector;

template IVector intervalHull<IVector>(const IVector &A, const IVector &B);
template void split<IVector,IVector>(IVector& v, IVector& rv);
template void split<IVector,DVector>(const IVector& v, DVector&, IVector& rv);
template IVector midVector<IVector>(const IVector& v);
template IVector leftVector<IVector>(const IVector& v);
template IVector rightVector<IVector>(const IVector& v);
template capd::DInterval maxDiam<IVector>(const IVector& v);
template double maxWidth<IVector>(const IVector& v);
template IVector diam<IVector>(const IVector& v);
template DVector widths<IVector>(const IVector& v);
template IVector intervalBall<IVector>(const IVector &iv, const capd::DInterval &r);
template capd::DInterval solveAffineInclusion<IVector>(const IVector& a,const IVector& p,const IVector& c);
template capd::DInterval solveAffineInclusion<IVector>(const IVector& a,const IVector& p,const IVector& c,int&);
template bool subset<IVector>(const IVector& v1, const IVector& v2);
template bool subsetInterior<IVector>(const IVector& v1, const IVector& v2);
template IVector intersection<IVector>(const IVector &v1, const IVector &v2);
template std::string vectorToString<IVector>( const IVector & v, int firstIndex , int lastIndex , int precision);
template std::ostream & printVector<IVector>(std::ostream & str, const IVector & v, int firstIndex, int lastIndex);
template std::ostream & printVector<>(std::ostream & str, const ColumnVector<capd::DInterval, CAPD_DEFAULT_DIMENSION> &, int, int);

template IVector multiplyObjectScalar<>(ColumnVector<capd::DInterval, CAPD_DEFAULT_DIMENSION>  const&, capd::DInterval const&);
template IVector multiplyObjectScalar<>(RowVector<capd::DInterval, CAPD_DEFAULT_DIMENSION>  const&, capd::DInterval const&);
template IVector multiplyObjectScalar<>(IVector  const&, capd::DInterval const&);
template IVector multiplyObjectScalar<>(IVector  const&, long const&);

template IVector divideObjectScalar<>(IVector const&, long const&);
template IVector divideObjectScalar<>(IVector const&, capd::DInterval const&);


template ColumnVector<capd::DInterval, CAPD_DEFAULT_DIMENSION>& divideAssignObjectScalar<>(ColumnVector<capd::DInterval, CAPD_DEFAULT_DIMENSION>&, capd::DInterval const&);
template RowVector<capd::DInterval, CAPD_DEFAULT_DIMENSION>& divideAssignObjectScalar<>(RowVector<capd::DInterval, CAPD_DEFAULT_DIMENSION>&, capd::DInterval const&);
template IVector& multiplyAssignObjectScalarAddObject<>(IVector&, capd::DInterval const&, IVector const&);
template void clear<>(ColumnVector<capd::DInterval, CAPD_DEFAULT_DIMENSION>&);
template bool equal<>(Vector<capd::DInterval, CAPD_DEFAULT_DIMENSION> const&, Vector<capd::DInterval, CAPD_DEFAULT_DIMENSION> const&);

 template Vector<capd::DInterval, CAPD_DEFAULT_DIMENSION>::ScalarType scalarProduct<>(Vector<capd::DInterval, CAPD_DEFAULT_DIMENSION> const&, Vector<capd::DInterval, CAPD_DEFAULT_DIMENSION> const&);

  }}  // end of namespace capd::vectalg
