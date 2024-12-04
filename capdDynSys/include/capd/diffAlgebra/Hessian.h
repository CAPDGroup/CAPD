/// @addtogroup diffAlgebra
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file Hessian.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DIFFALGEBRA_HESSIAN_H_
#define _CAPD_DIFFALGEBRA_HESSIAN_H_

#include "capd/vectalg/Container.h"
#include "capd/vectalg/Matrix.h"

namespace capd{
namespace diffAlgebra{

using capd::vectalg::__size_type;
using capd::vectalg::__difference_type;

/**
 * This class is used to store second order partial derivatives of a function \f$ R^D\to R^M \f$.
 *
 * It provides indexing, iterators and some algorithms for hessians.
 *
 */
template<typename Scalar, __size_type M, __size_type D>
class Hessian : public capd::vectalg::Container<Scalar,M*D*(1+D)/2>
{
public:
  typedef Scalar ScalarType;
  typedef capd::vectalg::Container<ScalarType,M*D*(1+D)/2> ContainerType;
  typedef typename ContainerType::iterator iterator;
  typedef typename ContainerType::const_iterator const_iterator;
  typedef __size_type size_type;
  typedef __difference_type difference_type;
  typedef capd::vectalg::ColumnVector<Scalar,M> RefColumnVectorType;

  // constructors
  Hessian();
  explicit Hessian(size_type dim);
  Hessian(size_type image, size_type domain);
  Hessian(size_type image, size_type domain, const ScalarType* data);

  ScalarType& operator()(size_type i, size_type j, size_type k); ///< value of\f$ \partial^2f_i/{\partial j\partial k}\f$
  const ScalarType& operator()(size_type i, size_type j, size_type k) const; ///< value of \f$\partial^2f_i/{\partial j\partial k}\f$
  RefColumnVectorType operator()(size_type j, size_type k); ///< returns vector of derivatives \f$\partial^2f/{\partial j\partial k}\f$

  //Hessian& operator=(const Hessian& h);     ///< assignment
  Hessian& operator+=(const Hessian& s);    ///< component-wise addition
  Hessian& operator*=(const ScalarType& s); ///< multiplication of each coefficient by a scalar

  size_type dimension() const;        ///< dimension of the phase space
  size_type imageDimension() const;   ///< dimension of codomain
  void resize(size_type _M, size_type _D) {
    ContainerType::resize(_M*_D*(1+_D)/2);
    this->m_imageDimension = M>0 ? M : _M;
    this->m_domainDimension = D>0 ? D : _D;
  }
  static Hessian* makeArray(size_type N, size_type _M, size_type _D);

  using ContainerType::operator[];
  using ContainerType::begin;
  using ContainerType::end;
  using ContainerType::clear;
protected:
  size_type m_domainDimension;
  size_type m_imageDimension;
}; // the end of class Hessian

template<typename Scalar,__size_type M, __size_type D>
inline bool operator== (const Hessian<Scalar,M,D>& h1, const Hessian<Scalar,M,D>& h2){
  return capd::vectalg::equal(h1,h2);
}


template<typename Scalar,__size_type M, __size_type D>
inline bool operator!= (const Hessian<Scalar,M,D>& h1, const Hessian<Scalar,M,D>& h2){
  return capd::vectalg::notEqual(h1,h2);
}


template<typename ScalarType, __size_type D>
Hessian<ScalarType,D,D> operator*(const capd::vectalg::Matrix<ScalarType,D,D>& m, const Hessian<ScalarType,D,D>& c2);

template<typename ScalarType, __size_type D>
Hessian<ScalarType,D,D> operator*(const Hessian<ScalarType,D,D>& c2, const capd::vectalg::Matrix<ScalarType,D,D>& m);

template<typename ScalarType, __size_type D>
Hessian<ScalarType,D,D> operator*(ScalarType c, const Hessian<ScalarType,D,D>& H);

template <typename ScalarType,__size_type D>
inline
Hessian<ScalarType,D,D> operator+(const Hessian<ScalarType,D,D>& H1, const Hessian<ScalarType,D,D>& H2){
  Hessian<ScalarType,D,D> result(H1.dimension());
  capd::vectalg::addObjects(H1,H2,result);
  return result;
}

// ------------------- inline definitions -------------------

/*
template<typename ScalarType, __size_type M, __size_type D>
inline Hessian<ScalarType,M,D>& Hessian<ScalarType,M,D>::operator=(const Hessian& h)
{
  ContainerType::operator=(h);
  return *this;
}
 */

template<typename ScalarType, __size_type M, __size_type D>
inline Hessian<ScalarType,M,D>::Hessian() : m_domainDimension(D), m_imageDimension(M)
{}

template<typename ScalarType, __size_type M, __size_type D>
inline Hessian<ScalarType,M,D>::Hessian(size_type dim)
  : ContainerType(dim*dim*(1+dim)/2)
{
  m_domainDimension = D>0 ? D : dim;
  m_imageDimension = M>0 ? M : dim;
}

template<typename ScalarType, __size_type M, __size_type D>
inline Hessian<ScalarType,M,D>::Hessian(size_type image, size_type domain)
  : ContainerType(image*domain*(1+domain)/2)
{
  m_domainDimension = D>0 ? D : domain;
  m_imageDimension = M>0 ? M : image;
}

template<typename ScalarType, __size_type M, __size_type D>
inline Hessian<ScalarType,M,D>::Hessian(size_type image, size_type domain, const ScalarType* data)
  : ContainerType(image*domain*(1+domain)/2,true)
{
  std::copy(data,data+this->size(),this->begin());
  m_domainDimension = D>0 ? D : domain;
  m_imageDimension = M>0 ? M : image;
}

template<typename ScalarType, __size_type M, __size_type D>
inline typename Hessian<ScalarType,M,D>::size_type Hessian<ScalarType,M,D>::dimension() const
{
  return m_domainDimension;
}

template<typename ScalarType, __size_type M, __size_type D>
inline typename Hessian<ScalarType,M,D>::size_type Hessian<ScalarType,M,D>::imageDimension() const
{
  return m_imageDimension;
}

template<typename ScalarType, __size_type M, __size_type D>
inline ScalarType& Hessian<ScalarType,M,D>::operator()(size_type i, size_type j, size_type k)
{
  return k<=j ?
      ContainerType::operator[](i*m_domainDimension*(m_domainDimension+1)/2 + j*(j+1)/2 + k)
    : ContainerType::operator[](i*m_domainDimension*(m_domainDimension+1)/2 + k*(k+1)/2 + j);
}

template<typename ScalarType, __size_type M, __size_type D>
inline const ScalarType& Hessian<ScalarType,M,D>::operator()(size_type i, size_type j, size_type k) const
{
  return k<=j ?
      ContainerType::operator[](i*m_domainDimension*(m_domainDimension+1)/2 + j*(j+1)/2 + k)
    : ContainerType::operator[](i*m_domainDimension*(m_domainDimension+1)/2 + k*(k+1)/2 + j);
}

template<typename ScalarType, __size_type M, __size_type D>
inline typename Hessian<ScalarType,M,D>::RefColumnVectorType
Hessian<ScalarType,M,D>::operator()(size_type j, size_type k)
{
  size_type pos = k<=j ? j*(j+1)/2 + k : k*(k+1)/2 + j;
  return RefColumnVectorType(this->begin()+pos,m_domainDimension*(m_domainDimension+1)/2,this->m_imageDimension);
}

template <typename T, capd::vectalg::__size_type M, capd::vectalg::__size_type D>
std::ostream & operator << (std::ostream & str, const Hessian<T,M,D> & h){
  str << "{";
  for(__size_type i=0; i<h.dimension(); i++){
    str << "{\n";
    for(__size_type j=0 ; j< h.imageDimension(); ++j){
      for(__size_type k=0; k <= j; ++k){
        str << h(i,j,k) <<" ";
      }
      str << "\n";
    }
    str<<(((i+1) != h.dimension())? "}," : "}");
  }
  str << "}" <<std::endl;
  return str;

}

}} // namespace capd::diffAlgebra

#endif // _CAPD_DIFFALGEBRA_C2COEFF_H_

/// @}
