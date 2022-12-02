
/////////////////////////////////////////////////////////////////////////////
/// @file vectalg/LDVector.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include <iostream>
#include "capd/vectalg/Dimension.h"
#include "capd/vectalg/Vector.hpp"

namespace capd{
  namespace vectalg{

template class Vector<long double,CAPD_DEFAULT_DIMENSION>;

template Vector<long double,CAPD_DEFAULT_DIMENSION> abs<long double,CAPD_DEFAULT_DIMENSION> (const Vector<long double,CAPD_DEFAULT_DIMENSION> &v);
template Vector<long double,CAPD_DEFAULT_DIMENSION> operator- <long double,CAPD_DEFAULT_DIMENSION> (const Vector<long double,CAPD_DEFAULT_DIMENSION> &v);
template Vector<long double,CAPD_DEFAULT_DIMENSION> operator+ <long double,CAPD_DEFAULT_DIMENSION> (const Vector<long double,CAPD_DEFAULT_DIMENSION> &v1,const Vector<long double,CAPD_DEFAULT_DIMENSION> &v2);
template Vector<long double,CAPD_DEFAULT_DIMENSION> operator- <long double,CAPD_DEFAULT_DIMENSION> (const Vector<long double,CAPD_DEFAULT_DIMENSION> &v1,const Vector<long double,CAPD_DEFAULT_DIMENSION> &v2);
template long double operator* <long double,CAPD_DEFAULT_DIMENSION> (const Vector<long double,CAPD_DEFAULT_DIMENSION> &v1,const Vector<long double,CAPD_DEFAULT_DIMENSION> &v2);

template Vector<long double,CAPD_DEFAULT_DIMENSION> operator* <long double,long double,CAPD_DEFAULT_DIMENSION> (const Vector<long double,CAPD_DEFAULT_DIMENSION> &v,const long double &s);
template Vector<long double,CAPD_DEFAULT_DIMENSION> operator* <long double,long double,CAPD_DEFAULT_DIMENSION> (const long double &s,const Vector<long double,CAPD_DEFAULT_DIMENSION> &v);
template Vector<long double,CAPD_DEFAULT_DIMENSION> operator/ <long double,long double,CAPD_DEFAULT_DIMENSION> (const Vector<long double,CAPD_DEFAULT_DIMENSION> &v, const long double &s);

template Vector<long double,CAPD_DEFAULT_DIMENSION> operator* <long double,long,CAPD_DEFAULT_DIMENSION>(const Vector<long double,CAPD_DEFAULT_DIMENSION> &v,const long &s);
template Vector<long double,CAPD_DEFAULT_DIMENSION> operator* <long double,long,CAPD_DEFAULT_DIMENSION>(const long &s,const Vector<long double,CAPD_DEFAULT_DIMENSION> &v);
template Vector<long double,CAPD_DEFAULT_DIMENSION> operator/ <long double,long,CAPD_DEFAULT_DIMENSION>(const Vector<long double,CAPD_DEFAULT_DIMENSION> &v, const long &s);

template Vector<long double,CAPD_DEFAULT_DIMENSION> operator* <long double,int,CAPD_DEFAULT_DIMENSION>(const Vector<long double,CAPD_DEFAULT_DIMENSION> &v,const int &s);
template Vector<long double,CAPD_DEFAULT_DIMENSION> operator* <long double,int,CAPD_DEFAULT_DIMENSION>(const int &s,const Vector<long double,CAPD_DEFAULT_DIMENSION> &v);
template Vector<long double,CAPD_DEFAULT_DIMENSION> operator/ <long double,int,CAPD_DEFAULT_DIMENSION>(const Vector<long double,CAPD_DEFAULT_DIMENSION> &v, const int &s);

template Vector<long double,CAPD_DEFAULT_DIMENSION> operator+ <long double,CAPD_DEFAULT_DIMENSION> (const Vector<long double,CAPD_DEFAULT_DIMENSION> &v,const long double &s);
template Vector<long double,CAPD_DEFAULT_DIMENSION> operator- <long double,CAPD_DEFAULT_DIMENSION> (const Vector<long double,CAPD_DEFAULT_DIMENSION> &v,const long double &s);
template bool operator < <long double,CAPD_DEFAULT_DIMENSION> (const Vector<long double,CAPD_DEFAULT_DIMENSION> &v1,const Vector<long double,CAPD_DEFAULT_DIMENSION> &v2);
template bool operator > <long double,CAPD_DEFAULT_DIMENSION> (const Vector<long double,CAPD_DEFAULT_DIMENSION> &v1,const Vector<long double,CAPD_DEFAULT_DIMENSION> &v2);
template bool operator<= <long double,CAPD_DEFAULT_DIMENSION> (const Vector<long double,CAPD_DEFAULT_DIMENSION> &v1,const Vector<long double,CAPD_DEFAULT_DIMENSION> &v2);
template bool operator>= <long double,CAPD_DEFAULT_DIMENSION> (const Vector<long double,CAPD_DEFAULT_DIMENSION> &v1,const Vector<long double,CAPD_DEFAULT_DIMENSION> &v2);
template std::ostream &operator<< <long double,CAPD_DEFAULT_DIMENSION> (std::ostream &out, const Vector<long double,CAPD_DEFAULT_DIMENSION> &v);
template std::istream &operator>> <long double,CAPD_DEFAULT_DIMENSION>(std::istream &inp, Vector<long double,CAPD_DEFAULT_DIMENSION> &v);


template void subtractObjects<>(const Vector<long double,CAPD_DEFAULT_DIMENSION>& v1,const Vector<long double,CAPD_DEFAULT_DIMENSION>& v2, Vector<long double,CAPD_DEFAULT_DIMENSION>& result);
template void addObjects<>(const Vector<long double,CAPD_DEFAULT_DIMENSION>& v1,const Vector<long double,CAPD_DEFAULT_DIMENSION>& v2, Vector<long double,CAPD_DEFAULT_DIMENSION>& result);
typedef Vector<long double,CAPD_DEFAULT_DIMENSION> LDVector;
template LDVector divideObjectScalar<>(const LDVector & v1, const long double & s);
//capd::vectalg::Vector<long double, 0u> capd::vectalg::divideObjectScalar<capd::vectalg::Vector<long double, 0u>, capd::vectalg::Vector<long double, 0u>, long double>(capd::vectalg::Vector<long double, 0u> const&, long double const&);
}}  // end of namespace capd::vectalg


