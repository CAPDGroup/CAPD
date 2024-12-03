/////////////////////////////////////////////////////////////////////////////
/// @file Vector_inline.h
///
/// This file contains inline definitions for class Vector
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.


#ifndef _CAPD_VECTALG_VECTOR_INLINE_H_ 
#define _CAPD_VECTALG_VECTOR_INLINE_H_ 

#include "capd/vectalg/Vector.h"
#include "capd/vectalg/algebraicOperations.h"

namespace capd{
namespace vectalg{

template<typename Scalar, __size_type dim>
class Vector;

template<typename Scalar, __size_type dim>
inline Vector<Scalar,dim> abs(const Vector<Scalar,dim>& v){
  return absoluteValue< Vector<Scalar,dim>, Vector<Scalar,dim> > (v);
}

template<typename Scalar, __size_type dim>
inline Vector<Scalar,dim> operator-(const Vector<Scalar,dim>& v){
  return unaryMinus< Vector<Scalar,dim>, Vector<Scalar,dim> >(v);
}

template<typename Scalar, __size_type dim>
inline Vector<Scalar,dim> operator+(const Vector<Scalar,dim>& v1,const Vector<Scalar,dim>& v2){
  return addObjects< Vector<Scalar,dim>, Vector<Scalar,dim>, Vector<Scalar,dim> >(v1,v2);  
}

template<typename Scalar, __size_type dim>
inline Vector<Scalar,dim> operator-(const Vector<Scalar,dim>& v1,const Vector<Scalar,dim>& v2){
  return subtractObjects< Vector<Scalar,dim>, Vector<Scalar,dim>, Vector<Scalar,dim> >(v1,v2);  
}

template<typename Scalar, __size_type dim>
inline Scalar operator*(const Vector<Scalar,dim>& v1,const Vector<Scalar,dim>& v2){
  return scalarProduct< Vector<Scalar,dim>, Vector<Scalar,dim> >(v1,v2);
}

template<typename Scalar, __size_type dim>
inline Vector<Scalar,dim> operator+(const Vector<Scalar,dim>& v, const Scalar& s){
  return addObjectScalar< Vector<Scalar,dim>, Vector<Scalar,dim>, Scalar >(v,s);
}

template<typename Scalar, __size_type dim>
inline Vector<Scalar,dim> operator-(const Vector<Scalar,dim>& v,const Scalar& s){
  return subtractObjectScalar< Vector<Scalar,dim>, Vector<Scalar,dim>, Scalar >(v,s);
}


template<typename Scalar, typename FactorType, __size_type dim>
inline Vector<Scalar,dim> operator*(const Vector<Scalar,dim>& v, const FactorType& s){
  return multiplyObjectScalar< Vector<Scalar,dim>, Vector<Scalar,dim>, FactorType >(v,s);
}

template<typename Scalar, typename FactorType, __size_type dim>
inline Vector<Scalar,dim> operator*(const FactorType& s,const Vector<Scalar,dim>& v){
  return multiplyObjectScalar< Vector<Scalar,dim>, Vector<Scalar,dim>, FactorType >(v,s);
}

template<typename Scalar, typename FactorType, __size_type dim>
Vector<Scalar,dim> operator/(const Vector<Scalar,dim>& v, const FactorType& s){
  return divideObjectScalar< Vector<Scalar,dim>, Vector<Scalar,dim>, FactorType >(v,s);
}


template<typename Scalar,__size_type dim>
inline Vector<Scalar,dim>& Vector<Scalar,dim>::operator=(const Vector<Scalar,dim>& v){
  ContainerType::operator=(v);
  return *this;
}

template<typename Scalar,__size_type dim>
inline Vector<Scalar,dim>::Vector(size_type A_dimension) : ContainerType(A_dimension)
{}

template<typename Scalar,__size_type dim>
inline Vector<Scalar,dim>::Vector(const Vector& A_vect) : ContainerType(A_vect)
{}

template<typename Scalar,__size_type dim>
inline Vector<Scalar,dim>::Vector(const Scalar& x,const Scalar& y,const Scalar& z) : ContainerType(3,true){
  (*this)[0]=x; (*this)[1]=y; (*this)[2]=z;
}

template<typename Scalar,__size_type dim>
inline Vector<Scalar,dim>::Vector(size_type a_Dimension,bool) : ContainerType(a_Dimension,true)
{}

//template<typename Scalar,__size_type dim>
//inline Vector<Scalar,dim>::Vector(Vector&& v) : ContainerType(std::move(v)) {
//}
//template<typename Scalar,__size_type dim>
//inline Vector<Scalar,dim> & Vector<Scalar,dim>::operator=(Vector && v) {
//  ContainerType::operator=(std::move(v));
//  return *this;
//}

template<typename Scalar,__size_type dim>
inline Vector<Scalar,dim>::Vector(std::initializer_list<ScalarType> l) 
  : ContainerType(l.size(), false) {
    if(l.size() == this->size())
  		std::copy(l.begin(), l.end(), this->begin());
    else
      throw std::range_error("Constructor of Vector with static size got initializer list with wrong size.");
}

template<typename Scalar,__size_type dim>
inline Vector<Scalar,dim>& Vector<Scalar,dim>::operator+=(const Vector<Scalar,dim>& v){
   return addAssignObjectObject < Vector<Scalar,dim>, Vector<Scalar,dim> > (*this,v);
}


template<typename Scalar,__size_type dim>
inline Vector<Scalar,dim>& Vector<Scalar,dim>::operator-=(const Vector& v){
   return subtractAssignObjectObject < Vector<Scalar,dim>, Vector<Scalar,dim> > (*this,v);
}

template<typename Scalar,__size_type dim>
inline Vector<Scalar,dim>& Vector<Scalar,dim>::operator=(const Scalar& s){
   return assignFromScalar< Vector<Scalar,dim>, Scalar > (*this,s);
}


template<typename Scalar,__size_type dim>
inline Vector<Scalar,dim>& Vector<Scalar,dim>::operator+=(const Scalar& s){
  return addAssignObjectScalar< Vector<Scalar,dim>, Scalar > (*this,s);
}

template<typename Scalar,__size_type dim>
inline Vector<Scalar,dim>& Vector<Scalar,dim>::operator-=(const Scalar& s){
  return subtractAssignObjectScalar< Vector<Scalar,dim>, Scalar > (*this,s);
}

template<typename Scalar,__size_type dim>
inline Vector<Scalar,dim>& Vector<Scalar,dim>::operator*=(const Scalar& s){
  return multiplyAssignObjectScalar< Vector<Scalar,dim>, Scalar > (*this,s);
}

template<typename Scalar,__size_type dim>
inline Vector<Scalar,dim>& Vector<Scalar,dim>::operator/=(const Scalar& s){
  return divideAssignObjectScalar< Vector<Scalar,dim>, Scalar > (*this,s);
}


template<typename Scalar,__size_type dim>
inline bool operator< (const Vector<Scalar,dim>& v1, const Vector<Scalar,dim>& v2){
  return lessThan(v1,v2);
}


template<typename Scalar,__size_type dim>
inline bool operator> (const Vector<Scalar,dim>& v1, const Vector<Scalar,dim>& v2){
  return greaterThan(v1,v2);
}


template<typename Scalar,__size_type dim>
inline bool operator<= (const Vector<Scalar,dim>& v1, const Vector<Scalar,dim>& v2){
  return lessEqual(v1,v2);
}


template<typename Scalar,__size_type dim>
inline bool operator>= (const Vector<Scalar,dim>& v1, const Vector<Scalar,dim>& v2){
  return greaterEqual(v1,v2);
}

// ----------------------------------- equality --------------------------------------

template<typename Scalar,__size_type dim>
inline bool operator== (const Vector<Scalar,dim>& v1, const Vector<Scalar,dim>& v2){
  return equal(v1,v2);
}


template<typename Scalar,__size_type dim>
inline bool operator!= (const Vector<Scalar,dim>& v1, const Vector<Scalar,dim>& v2){
  return notEqual(v1,v2);
} 

template<typename Scalar,__size_type dim>
inline typename Vector<Scalar,dim>::ScalarType Vector<Scalar,dim>::euclNorm(void) const{
  return capd::vectalg::euclNorm(*this);
}

template<typename Scalar,__size_type dim>
inline bool Vector<Scalar,dim>::normalize(){
  return capd::vectalg::normalize(*this);
}

template<typename Scalar,__size_type dim>
template<typename S,
        typename std::enable_if<std::is_convertible<S, Scalar>::value && !std::is_same<S, Scalar>::value , int>::type>
inline Vector<Scalar,dim>::Vector(const Vector<S,dim>& v) : ContainerType(v.dimension(),true){
  assignObjectObject(*this,v);
}

}} // namespace capd::vectalg

#endif // _CAPD_VECTALG_VECTOR_INLINE_H_ 

