/// @addtogroup diffAlgebra
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file CnContainer.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include <stdexcept>
#include "capd/basicalg/factrial.h"
#include "capd/vectalg/Multiindex.h"
#include "capd/vectalg/Container.h"

#ifndef _CAPD_DIFFALGEBRA_CNCONTAINER_H_
#define _CAPD_DIFFALGEBRA_CNCONTAINER_H_

namespace capd{
namespace diffAlgebra{

using capd::vectalg::__size_type;
using capd::vectalg::__difference_type;

/**
 * The class is used to store coefficients of a multivariate polynomial of degree D
 * \f$ f:R^N->R^M \f$
 * Coefficients themselves can be polynomials as well.
 * Usually Object = Scalar or an univariate polynomial
 *
 * The total number of coefficients is equal to \f$ M {N+D\choose D} \f$
*/

template<typename Object, __size_type M, __size_type N, __size_type D>
class CnContainer : public capd::vectalg::Container<Object,M*N*D!=0 ? M*Binomial<N+D,D>::value : 0>
{
public:
  typedef capd::vectalg::Container<Object,M*N*D!=0 ? M*Binomial<N+D,D>::value : 0> BaseContainer;
  typedef Object ObjectType;
  typedef Object* iterator;
  typedef const Object* const_iterator;
  typedef capd::vectalg::Multipointer Multipointer;
  typedef capd::vectalg::Multiindex Multiindex;
  typedef __size_type size_type;
  typedef __difference_type difference_type;

  CnContainer(size_type m, size_type n, size_type d, const Object& p); ///< creates a container for polynomial of n variables, m components and degree d. Each element will be set to p.
  CnContainer(size_type m, size_type n, size_type d);                  ///< creates a container for polynomial of n variables, m components and degree d. Default constructor will be used to initialize each element in the container.
  CnContainer& operator=(const Object& p);           ///< assigns object p to each element of the container

  CnContainer(const CnContainer& v) = default;
  CnContainer(CnContainer&& v) : BaseContainer(std::move(v)), m_N(v.m_N), m_M(v.m_M), m_D(v.m_D) {}
  CnContainer & operator=(const CnContainer & v) = default;
  CnContainer & operator=(CnContainer && v) {
    swap(*this, v);
    return *this;
  }

  size_type imageDimension() const; ///< returns number of polynomials (dimension of counterdomain)
  size_type dimension() const;      ///< returns number of variables of the polynomial
  size_type degree() const;         ///< returns degree of the polynomial

//  void resize(int newDegree, bool copyData = true);                   ///< changes degree of the polynomial
//  void resize(int newDegree, int newDimension, bool copyData = true); ///< changes degree and the number of variables of the polynomial
//  void resize(int newDegree, int newDimension, int newimageDim, bool copyData = true); ///< changes degree and the number of variables of the polynomial

// indexing
  
  using BaseContainer::operator[]; //< direct access to an element by its absolute position in the container

  Object& operator()(size_type i, const Multipointer& mp);  ///< selection of coefficient of i-th component that correspond to multipointer mp
  Object& operator()(size_type i, const Multipointer&, const Multipointer&);
  Object& operator()(size_type i, const Multiindex& mi); ///< selection of coefficient of i-th component that correspond to multiindex mi

  const Object& operator()(size_type i, const Multipointer&) const; ///< selection of coefficient of i-th component that correspond to multipointer mp
  const Object& operator()(size_type i, const Multipointer&, const Multipointer&) const;
  const Object& operator()(size_type i, const Multiindex&) const; ///< selection of coefficient of i-th component that correspond to multiindex mi

// operators for C^0, C^1, C^2 and C^3 algorithms
  Object& operator()(size_type i);                      ///< returns constant term of the i-th component of polynomial
  Object& operator()(size_type i, size_type j);               ///< returns reference to a coefficient in linear part, i.e. \f$ df_i/dx_j \f$
  Object& operator()(size_type i, size_type j, size_type c);        ///< returns reference to a coefficient in second order part, i.e. \f$ d^2f_i/dx_jdx_c \f$
  Object& operator()(size_type i, size_type j, size_type c, size_type k); ///< returns reference to a coefficient in third order part, i.e. \f$ d^3f_i/dx_jdx_cdx_k \f$

  const Object& operator()(size_type i) const;                      ///< returns constant term of the i-th component of polynomial
  const Object& operator()(size_type i, size_type j) const;               ///< returns read only reference to a coefficient in linear part, i.e. \f$ df_i/dx_j \f$
  const Object& operator()(size_type i, size_type j, size_type c) const;        ///< returns read only reference to a coefficient in second order part, i.e. \f$ d^2f_i/dx_jdx_c \f$
  const Object& operator()(size_type i, size_type j, size_type c, size_type k) const; ///< returns read only reference to a coefficient in third order part, i.e. \f$ d^3f_i/dx_jdx_cdx_k \f$
   
// iterators
  using BaseContainer::begin;   //< iterator selection. Returns iterator to the first element in container
  using BaseContainer::end;     //< iterator selection. Returns iterator to the first element in container
  using BaseContainer::clear;

  iterator begin(size_type i);        ///< iterator selection. Returns iterator to the first coefficient of the i-th component
  iterator end(size_type i);          ///< iterator selection. Returns iterator to an element after the last element the i-th component
  iterator begin(size_type i, size_type d); ///< iterator selection. Returns iterator to the first coefficient of the i-th component of the homogeneous part of degree 'd'
  iterator end(size_type i, size_type d);   ///< iterator selection. Returns iterator to an element after the last coefficient of the i-th component of the homogeneous part of degree 'd'

  const_iterator begin(size_type i) const;            ///< iterator selection. Returns iterator to the first coefficient of the i-th component
  const_iterator end(size_type i) const;              ///< iterator selection. Returns iterator to an element after the last element the i-th component
  const_iterator begin(size_type i, size_type d) const;     ///< iterator selection. Returns iterator to the first coefficient of the i-th component of the homogeneous part of degree 'd'
  const_iterator end(size_type i, size_type d) const;       ///< iterator selection. Returns iterator to an element after the last coefficient of the i-th component of the homogeneous part of degree 'd'

/**
 * Selection of elements by multipointers.
 *
 * Iterators do not give information about the index of partial derivative.
 * Access by multipointer is significantly slower than by iterator because the multipointer must be recomputed to the index in array.
 *
 * Typical usage of multipointers is as follows:
 * <code>
 * Multipointer mp = cnContainer.first(d);
 * int i = ...; // fix i-th component
 * do{
 *   // do something
 *   cout << mp << "\t" << cnContainer(i,mp) << endl;
 * }while(cnContainer.hasNext(mp));
 * </code>
 * Iterators and multipointers read coefficients of a homogeneous polynomial in the same order.
 */
  Multipointer first(size_type d) const;
  bool hasNext(Multipointer&) const; ///< see description of the method first.
  bool hasNext(Multiindex&) const; ///< see description of the method first.

  friend void swap(CnContainer & A, CnContainer & B){ ///< swaps the content of two containers
    std::swap(A.m_N, B.m_N);
    std::swap(A.m_M, B.m_M);
    std::swap(A.m_D, B.m_D);
    swap(static_cast<BaseContainer &>(A),static_cast<BaseContainer &>(B));
  }

protected:
  size_type m_N; ///< number of variables
  size_type m_M; ///< number of components
  size_type m_D; ///< total degree of polynomial
}; // the end of class CnContainer

// ------------------- member definitions -------------------
// the following function computes the next multipointer after mp
// it returns false if mp is the last multipointer on a given level

template<typename Object, __size_type M, __size_type N, __size_type D>
bool CnContainer<Object,M,N,D>::hasNext(Multipointer& mp) const
{
  return mp.hasNext(this->dimension());
}

template<typename Object, __size_type M, __size_type N, __size_type D>
bool CnContainer<Object,M,N,D>::hasNext(Multiindex& mp) const
{
  return mp.hasNext();
}

// ------------------- member definitions -------------------

template<typename Object, __size_type M, __size_type N, __size_type D>
CnContainer<Object,M,N,D>& CnContainer<Object,M,N,D>::operator=(const Object& p)
{
  iterator b=begin(), e=end();
  while(b!=e)
  {
    *b=p;
    ++b;
  }
  return *this;
}

// ----------------------------------------------------------

template<typename Object, __size_type M, __size_type N, __size_type D>
inline CnContainer<Object,M,N,D>::CnContainer(size_type m, size_type n, size_type d, const Object& p)
  : BaseContainer(m*binomial(n+d,d),true)
{
  m_N = N>0 ? N : n;
  m_M = M>0 ? M : m;
  m_D = D>0 ? D : d;

  iterator b = this->begin(), e=this->end();
  while(b!=e)
  {
    *b = p;
    ++b;
  }
}

// ----------------------------------------------------------
/*
template<typename Object, int M, int N, int D>
void CnContainer<Object,M,N,D>::resize(int newDegree,  bool copyData)
{
  if(newDegree == m_D)
    return;

  BaseContainer::resize(m_M*newton(m_N+m_D,m_D));
  if(copyData){
    int minDegree = m_D < newDegree ? m_D : newDegree;
    for(int i=0;i<m_N;++i)
    {
      iterator b = begin(i), e=end(i,minDegree);
      Object* p =  newData+ i*newton(m_dim+newRank,m_dim);
      while(b!=e)
      {
        *p = *b;
        ++p;
        ++b;
      }
    }
  }
  delete [] m_data;
  m_data = newData;
  m_rank = newRank;
  m_size = newSize;
}
*/
// ----------------------------------------------------------

/**
 * Resizes CnContainer
 * @param newRank       new maximal order
 * @param newDimension  new dimension
 * @param copyData      flag that controls if data is copied
 */
/*
template<typename Object, int M, int N, int D>
void CnContainer<Object,M,N,D>::resize(int newRank, int newDimension, bool copyData)
{
  if((newRank == m_rank) && (newDimension == m_dim))
    return;

  int newSize = newDimension*newton(newDimension+newRank,newDimension);

  Object* newData = new Object[newSize];
  if(copyData){
    int minRank = m_rank < newRank ? m_rank : newRank;
    int minDim = m_dim < newDimension ? m_dim : newDimension;
    for(int i=0; i< minDim; ++i){
      iterator b = begin(i),
          e = end(i,minRank);
      Object* p =  newData+ i*newton(newDimension+newRank,newDimension);
      while(b!=e)
      {
        *p = *b;
        ++p;
        ++b;
      }
    }
  }
  delete [] m_data;
  m_data = newData;
  m_rank = newRank;
  m_size = newSize;
  m_dim = newDimension;
}
*/
// ------------------- inline definitions -------------------

template<typename Object, __size_type M, __size_type N, __size_type D>
inline typename CnContainer<Object,M,N,D>::Multipointer CnContainer<Object,M,N,D>::first(size_type degree) const
{
  return Multipointer(degree);
}

// ----------------------------------------------------------

template<typename Object, __size_type M, __size_type N, __size_type D>
inline CnContainer<Object,M,N,D>::CnContainer(size_type m, size_type n, size_type d)
  : BaseContainer(m*binomial(n+d,d))
{
  m_N = N>0 ? N : n;
  m_M = M>0 ? M : m;
  m_D = D>0 ? D : d;
}

// ----------------------------------------------------------

template<typename Object, __size_type M, __size_type N, __size_type D>
inline typename CnContainer<Object,M,N,D>::size_type CnContainer<Object,M,N,D>::dimension() const
{
  return m_N;
}

// ----------------------------------------------------------

template<typename Object, __size_type M, __size_type N, __size_type D>
inline typename CnContainer<Object,M,N,D>::size_type CnContainer<Object,M,N,D>::imageDimension() const
{
  return m_M;
}

// ----------------------------------------------------------

template<typename Object, __size_type M, __size_type N, __size_type D>
inline typename CnContainer<Object,M,N,D>::size_type CnContainer<Object,M,N,D>::degree() const
{
  return m_D;
}

// ----------------------------------------------------------

template<typename Object, __size_type M, __size_type N, __size_type D>
inline Object& CnContainer<Object,M,N,D>::operator()(size_type i, const Multipointer& mp)
{
  return *(begin(i,mp.dimension()) + mp.index(m_N,m_D));
}

// ----------------------------------------------------------

template<typename Object, __size_type M, __size_type N, __size_type D>
inline const Object& CnContainer<Object,M,N,D>::operator()(size_type i, const Multipointer& mp) const
{
  return *(begin(i,mp.dimension()) + mp.index(m_N,m_D));
}

// ----------------------------------------------------------

template<typename Object, __size_type M, __size_type N, __size_type D>
inline Object& CnContainer<Object,M,N,D>::operator()(size_type i, const Multipointer& mp, const Multipointer& sub)
{
  return *(begin(i,sub.dimension()) + mp.index(m_N,m_D,sub));
}

// ----------------------------------------------------------

template<typename Object, __size_type M, __size_type N, __size_type D>
inline const Object& CnContainer<Object,M,N,D>::operator()(size_type i, const Multipointer& mp, const Multipointer& sub) const
{
  return *(begin(i,sub.dimension()) + mp.index(m_N,m_D,sub));
}

// ----------------------------------------------------------

template<typename Object, __size_type M, __size_type N, __size_type D>
inline Object& CnContainer<Object,M,N,D>::operator()(size_type i, const Multiindex& mi)
{
  if(this->dimension()!=mi.dimension())
    throw std::runtime_error("CnContainer::operator(int,Multiindex) - incompatible dimensions of CnContainer and Multiindex");
  return *(begin(i,mi.module()) + mi.index(m_D));
}

// ----------------------------------------------------------

template<typename Object, __size_type M, __size_type N, __size_type D>
inline const Object& CnContainer<Object,M,N,D>::operator()(size_type i, const Multiindex& mi) const
{
  if(this->dimension()!=mi.dimension())
    throw std::runtime_error("CnContainer::operator(int,Multiindex) - incompatible dimensions of CnContainer and Multiindex");
  return *(begin(i,mi.module()) + mi.index(m_D));
}

// -------------- indexing for C^0-C^3 algorithms ------------

template<typename Object, __size_type M, __size_type N, __size_type D>
inline Object& CnContainer<Object,M,N,D>::operator()(size_type i)
{
  return *(this->begin() + i*binomial(m_D+m_N,m_D));
}

// ----------------------------------------------------------

template<typename Object, __size_type M, __size_type N, __size_type D>
inline const Object& CnContainer<Object,M,N,D>::operator()(size_type i) const
{
  return *(this->begin() + i*binomial(m_D+m_N,m_D));
}

// ----------------------------------------------------------

template<typename Object, __size_type M, __size_type N, __size_type D>
inline Object& CnContainer<Object,M,N,D>::operator()(size_type i, size_type j)
{
  return *(begin(i,1)+j);
}

// ----------------------------------------------------------

template<typename Object, __size_type M, __size_type N, __size_type D>
inline const Object& CnContainer<Object,M,N,D>::operator()(size_type i, size_type j) const
{
  return *(begin(i,1)+j);
}

// ----------------------------------------------------------

template<typename Object, __size_type M, __size_type N, __size_type D>
inline Object& CnContainer<Object,M,N,D>::operator()(size_type i, size_type j, size_type c)
{
  return j<=c ? 
    *(begin(i,2)+c-(j*(j+1))/2+j*m_N) :
    *(begin(i,2)+j-(c*(c+1))/2+c*m_N);
}

// ----------------------------------------------------------

template<typename Object, __size_type M, __size_type N, __size_type D>
inline const Object& CnContainer<Object,M,N,D>::operator()(size_type i, size_type j, size_type c) const
{
  return j<=c ? 
    *(begin(i,2)+c-(j*(j+1))/2+j*m_N) :
    *(begin(i,2)+j-(c*(c+1))/2+c*m_N);
}

// ----------------------------------------------------------

template<typename Object, __size_type M, __size_type N, __size_type D>
inline Object& CnContainer<Object,M,N,D>::operator()(size_type i, size_type j, size_type c, size_type k)
{
  // assume j<=c<=k
  if(c<j || k<c)
    throw std::runtime_error("CnContainer::operator(int,int,int,int) - incorrect indexes");
  return *(
    begin(i,3) +
    (j*( (j-1)*(j-2) + 3*m_N*(m_N-j+2) ))/6    +   ((j-c)*(c+j-2*m_N-1))/2 + k-c
  );
}

// ----------------------------------------------------------

template<typename Object, __size_type M, __size_type N, __size_type D>
inline const Object& CnContainer<Object,M,N,D>::operator()(size_type i, size_type j, size_type c, size_type k) const
{
  // assume j<=c<=k
  if(c<j || k<c)
    throw std::runtime_error("CnContainer::operator(int,int,int,int) - incorrect indexes");
  return *(
    begin(i,3) +
    (j*( (j-1)*(j-2) + 3*m_N*(m_N-j+2) ))/6    +   ((j-c)*(c+j-2*m_N-1))/2 + k-c
  );
}

// ----------------------------------------------------------

template<typename Object, __size_type M, __size_type N, __size_type D>
inline typename CnContainer<Object,M,N,D>::iterator CnContainer<Object,M,N,D>::begin(size_type i)
{
  return iterator(this->begin()+i*binomial(m_N+m_D,m_D));
}

// ----------------------------------------------------------

template<typename Object, __size_type M, __size_type N, __size_type D>
inline typename CnContainer<Object,M,N,D>::iterator CnContainer<Object,M,N,D>::end(size_type i)
{
  return iterator(this->begin()+(i+1)*binomial(m_N+m_D,m_D));
}

// ----------------------------------------------------------

template<typename Object, __size_type M, __size_type N, __size_type D>
inline typename CnContainer<Object,M,N,D>::iterator CnContainer<Object,M,N,D>::begin(size_type i, size_type degree)
{
  return degree==0
         ? iterator(this->begin()+ i*binomial(m_N+m_D,m_D))
         : iterator(this->begin()+ i*binomial(m_N+m_D,m_D) + binomial(m_N+degree-1,m_N));
}

// ----------------------------------------------------------

template<typename Object, __size_type M, __size_type N, __size_type D>
inline typename CnContainer<Object,M,N,D>::iterator CnContainer<Object,M,N,D>::end(size_type i, size_type degree)
{
  return iterator(this->begin() + i*binomial(m_N+m_D,m_D) + binomial(m_N+degree,m_N));
}

// ----------------------------------------------------------

template<typename Object, __size_type M, __size_type N, __size_type D>
inline typename CnContainer<Object,M,N,D>::const_iterator CnContainer<Object,M,N,D>::begin(size_type i) const
{
  return const_iterator(this->begin()+i*binomial(m_N+m_D,m_N));
}

// ----------------------------------------------------------

template<typename Object, __size_type M, __size_type N, __size_type D>
inline typename CnContainer<Object,M,N,D>::const_iterator CnContainer<Object,M,N,D>::end(size_type i) const
{
  return const_iterator(this->begin()+(i+1)*binomial(m_N+m_D,m_D));
}

// ----------------------------------------------------------

template<typename Object, __size_type M, __size_type N, __size_type D>
inline typename CnContainer<Object,M,N,D>::const_iterator CnContainer<Object,M,N,D>::begin(size_type i, size_type degree) const
{
  return degree==0
         ? const_iterator(this->begin() + i*binomial(m_N+m_D,m_N))
         : const_iterator(this->begin() + i*binomial(m_N+m_D,m_N) + binomial(m_N+degree-1,m_N));
}

// ----------------------------------------------------------

template<typename Object, __size_type M, __size_type N, __size_type D>
inline typename CnContainer<Object,M,N,D>::const_iterator CnContainer<Object,M,N,D>::end(size_type i, size_type degree) const
{
  return const_iterator(this->begin() + i*binomial(m_N+m_D,m_N) + binomial(m_N+degree,m_N));
}

/// checks if two CnContainers are exactly the same.
template<typename Object, __size_type M, __size_type N, __size_type D>
bool operator == (const CnContainer<Object,M,N,D> & c1, const CnContainer<Object,M,N,D> & c2 ){
  if((c1.degree()()!=c2.degree()()) || (c1.dimension()!=c2.dimension()))
    return false;
  typename CnContainer<Object,M,N,D>::const_iterator it_c1 = c1.begin(), it_c2 = c2.begin(), end_c1 = c1.end();
  while(it_c1!=end_c1){
    if(*it_c1 != *it_c2)
       return false;
    it_c1++; it_c2++;
  }
  return true;
}

}} //namespace capd::diffAlgebra

#endif // _CAPD_DIFFALGEBRA_CNCONTAINER_H_

/// @}
