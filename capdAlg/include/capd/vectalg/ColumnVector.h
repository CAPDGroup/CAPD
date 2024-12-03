/// @addtogroup vectalg
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file ColumnVector.h
///
/// This file provides a template class ColumnVector
/// This class realizes a vector without own container, which is a reference
/// to a subset of other object with his own container.
/// A typical situation is a column of matrix which can be consider as a vector
///
/// The file 'RowVector.h' defines similar class, but in that case it is assumed
/// that data fill continuous part of a memory
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_VECTALG_COLUMNVECTOR_H_
#define _CAPD_VECTALG_COLUMNVECTOR_H_

#include <ostream>
#include "capd/vectalg/Vector.h"

namespace capd{
namespace vectalg{

template<typename Scalar>
class ColumnIterator{
public:
  typedef Scalar ScalarType;
  typedef __difference_type difference_type;

  inline ColumnIterator(ScalarType *p, difference_type pointerStride)
  : m_Pointer(p), m_pointerStride(pointerStride)
  {}

  inline ColumnIterator operator++(int)
  {
    ColumnIterator it(*this);
    m_Pointer += m_pointerStride;
    return it;
  }

  inline ColumnIterator operator--(int)
  {
    ColumnIterator it(*this);
    m_Pointer -= m_pointerStride;
    return it;
  }

  inline ColumnIterator& operator++()
  {
    m_Pointer += m_pointerStride;
    return *this;
  }

  inline ColumnIterator& operator--()
  {
    m_Pointer -= m_pointerStride;
    return *this;
  }

  inline ColumnIterator& operator+=(difference_type jump)
  {
    m_Pointer+= jump*m_pointerStride;
    return *this;
  }

  inline ColumnIterator& operator-=(difference_type jump)
  {
    m_Pointer-= jump*m_pointerStride;
    return *this;
  }

  inline bool operator!=(const ColumnIterator& second)
  {
    return m_Pointer!=second.m_Pointer;
  }

  inline ScalarType& operator*()
  {
    return *m_Pointer;
  }

  inline ScalarType* operator->()
  {
    return m_Pointer;
  }

private:
  ScalarType *m_Pointer;
  difference_type m_pointerStride;
  ColumnIterator(){} // we do not need a default constructor
};

// --------------------- const_iterator -------------------

template<typename Scalar>
class const_ColumnIterator{
public:
  typedef Scalar ScalarType;
  typedef __difference_type difference_type;

  inline const_ColumnIterator(const ScalarType* p, difference_type pointerStride)
    : m_Pointer(p), m_pointerStride(pointerStride)
  {}

  inline const_ColumnIterator operator++(int)
  {
    const_ColumnIterator it(*this);
    m_Pointer += m_pointerStride;
    return it;
  }

  inline const_ColumnIterator operator--(int)
  {
    const_ColumnIterator it(*this);
    m_Pointer -= m_pointerStride;
    return it;
  }

  inline const_ColumnIterator& operator++()
  {
    m_Pointer += m_pointerStride;
    return *this;
  }

  inline const_ColumnIterator& operator--()
  {
    m_Pointer -= m_pointerStride;
    return *this;
  }

  inline const_ColumnIterator& operator+=(difference_type jump)
  {
    m_Pointer += jump*m_pointerStride;
    return *this;
  }

  inline const_ColumnIterator& operator-=(difference_type jump)
  {
    m_Pointer -= jump*m_pointerStride;
    return *this;
  }

  inline bool operator!=(const const_ColumnIterator& second)
  {
    return m_Pointer!=second.m_Pointer;
  }

  inline const ScalarType& operator*()
  {
    return *m_Pointer;
  }

  inline const ScalarType* operator->()
  {
    return m_Pointer;
  }
private:
  const ScalarType*  m_Pointer;
  const_ColumnIterator(){} // we do not need a default constructor
  difference_type m_pointerStride;
};


template<typename Scalar, __size_type rows, __size_type cols>
class Matrix;

template<typename Scalar,__size_type dim>
class Vector;

/// This class realizes a vector without its own container, which is a reference
/// to a subset of other object with his own container.
/// A typical situation is a column of matrix which can be considered as a vector
///
template<typename Scalar, __size_type rows>
class ColumnVector
{
public:
  typedef Scalar ScalarType;
  typedef capd::vectalg::ColumnIterator<Scalar> iterator;
  typedef capd::vectalg::const_ColumnIterator<Scalar> const_iterator;
  typedef ColumnVector VectorType;
  typedef ColumnVector ContainerType;
  typedef __size_type size_type;
  typedef __difference_type difference_type;

  ColumnVector(const Scalar* pointer, difference_type stride, size_type dim);
  ColumnVector& operator=(const ColumnVector&);
  ColumnVector& operator=(const Vector<Scalar,rows>&);
  ColumnVector& operator+=(const ColumnVector&);
  ColumnVector& operator+=(const Vector<Scalar,rows>&);
  ColumnVector& operator-=(const ColumnVector&);
  ColumnVector& operator-=(const Vector<Scalar,rows>&);
  ColumnVector& operator*=(const Scalar&);
  ColumnVector& operator/=(const Scalar&);
  operator Vector<Scalar,rows>() const;

  Scalar& operator[](size_type row);
  const Scalar& operator[](size_type row) const;

  Scalar euclNorm() const;
  bool normalize();
  void clear();
  size_type dimension() const;
  void next();

  iterator begin();
  iterator end();
  const_iterator begin() const;
  const_iterator end() const;

  void assertEqualSize(const ColumnVector& c) const;
protected:
  Scalar* m_pointer;
  difference_type m_stride;
  size_type m_dim;
}; // end of class ColumnVector

// -------------------------- inline definitions ------------------------

template<typename Scalar, __size_type rows>
inline std::ostream& operator << (std::ostream& out, const ColumnVector<Scalar,rows>& s){
  return out << Vector<Scalar,rows>(s);
}

template<typename Scalar, __size_type rows>
inline Scalar ColumnVector<Scalar,rows>::euclNorm() const{
  return capd::vectalg::euclNorm(*this);
}

template<typename Scalar, __size_type rows>
inline bool ColumnVector<Scalar,rows>::normalize(){
  return capd::vectalg::normalize(*this);
}

template<typename Scalar, __size_type rows>
inline ColumnVector<Scalar,rows>& ColumnVector<Scalar,rows>::operator=(const ColumnVector& v){
  return assignObjectObject(*this,v);
}

template<typename Scalar, __size_type rows>
inline ColumnVector<Scalar,rows>& ColumnVector<Scalar,rows>::operator=(const Vector<Scalar,rows>& v){
  return assignObjectObject(*this,v);
}

// ----------------------- iterator selection ---------------------------

template<typename Scalar, __size_type rows>
inline typename ColumnVector<Scalar,rows>::size_type
ColumnVector<Scalar,rows>::dimension() const{
  return m_dim;
}

template<typename Scalar, __size_type rows>
inline typename ColumnVector<Scalar,rows>::iterator
ColumnVector<Scalar,rows>::begin(){
  return iterator(m_pointer,m_stride);
}

template<typename Scalar, __size_type rows>
inline typename ColumnVector<Scalar,rows>::iterator
ColumnVector<Scalar,rows>::end(){
  return iterator(m_pointer+m_dim*m_stride, m_stride);
}

template<typename Scalar, __size_type rows>
inline typename ColumnVector<Scalar,rows>::const_iterator
ColumnVector<Scalar,rows>::begin() const{
  return const_iterator(m_pointer, m_stride);
}

template<typename Scalar, __size_type rows>
inline typename ColumnVector<Scalar,rows>::const_iterator
ColumnVector<Scalar,rows>::end() const{
  return const_iterator(m_pointer+m_dim*m_stride, m_stride);
}

// ------------------------------ resize -----------------------------------

template<typename Scalar, __size_type rows>
inline void ColumnVector<Scalar,rows>::assertEqualSize(const ColumnVector& c) const{
  if(m_dim!=c.dimension())
    throw std::runtime_error("Unequal dimensions in ColumnVector::assertEqualSize");
}

// ------------------------------ constructor -----------------------------

template<typename Scalar, __size_type rows>
inline ColumnVector<Scalar,rows>::ColumnVector(const Scalar* pointer, difference_type stride, size_type dim)
    : m_pointer(const_cast<Scalar*>(pointer)),
      m_stride(stride), m_dim(dim)
{}

template<typename Scalar, __size_type rows>
inline Scalar& ColumnVector<Scalar,rows>::operator[](size_type row){
  return *(m_pointer + row*m_stride);
}

template<typename Scalar, __size_type rows>
inline const Scalar& ColumnVector<Scalar,rows>::operator[](size_type row) const{
  return *(m_pointer + row*m_stride);
}

template<typename Scalar, __size_type rows>
void ColumnVector<Scalar,rows>::next(){
  m_pointer++;
}

// -------------------- operator + ------------------------------------------

template<typename Scalar, __size_type rows>
inline Vector<Scalar,rows> operator+(const Vector<Scalar,rows>& u,
                                       const ColumnVector<Scalar,rows>& v
                                      )
{
  return addObjects< Vector<Scalar,rows>, Vector<Scalar,rows>, ColumnVector<Scalar,rows> > (u,v);
}

template<typename Scalar, __size_type rows>
inline Vector<Scalar,rows> operator+(const ColumnVector<Scalar,rows>&v,
                                       const Vector<Scalar,rows>&u
                                      )
{
  return addObjects< Vector<Scalar,rows>, Vector<Scalar,rows>, ColumnVector<Scalar,rows> > (u,v);
}

template<typename Scalar, __size_type rows>
inline Vector<Scalar,rows> operator+(const ColumnVector<Scalar,rows>&u,
                                       const ColumnVector<Scalar,rows>&v
                                      )
{
  return addObjects< Vector<Scalar,rows>, ColumnVector<Scalar,rows>, ColumnVector<Scalar,rows> > (u,v);
}

// -------------------- operator - ------------------------------------------

template<typename Scalar, __size_type rows>
inline Vector<Scalar,rows> operator-(const Vector<Scalar,rows>& u,
                                       const ColumnVector<Scalar,rows>& v
                                      )
{
  return subtractObjects< Vector<Scalar,rows>, Vector<Scalar,rows>, ColumnVector<Scalar,rows> > (u,v);
}

template<typename Scalar, __size_type rows>
inline Vector<Scalar,rows> operator-(const ColumnVector<Scalar,rows>&v,
                                       const Vector<Scalar,rows>&u
                                      )
{
  return subtractObjects< Vector<Scalar,rows>, Vector<Scalar,rows>, ColumnVector<Scalar,rows> > (u,v);
}

template<typename Scalar, __size_type rows>
inline Vector<Scalar,rows> operator-(const ColumnVector<Scalar,rows>&u,
                                       const ColumnVector<Scalar,rows>&v
                                      )
{
  return subtractObjects< Vector<Scalar,rows>, ColumnVector<Scalar,rows>, ColumnVector<Scalar,rows> > (u,v);
}

// ------------------------------- scalar product ----------------------------

template<typename Scalar, __size_type rows>
inline Scalar operator*(const Vector<Scalar,rows>& u, const ColumnVector<Scalar,rows>& v){
  return scalarProduct(u,v);
}

template<typename Scalar, __size_type rows>
inline Scalar operator*(const ColumnVector<Scalar,rows>&v, const Vector<Scalar,rows>&u){
  return scalarProduct(u,v);
}

template<typename Scalar, __size_type rows>
inline Scalar operator*(const ColumnVector<Scalar,rows>&v, const ColumnVector<Scalar,rows>&u){
  return scalarProduct(u,v);
}

// ------------------------- unary minus ----------------------------------------

template<typename Scalar, __size_type rows>
inline Vector<Scalar,rows> operator-(const ColumnVector<Scalar,rows>&u){
  return unaryMinus< Vector<Scalar,rows>, ColumnVector<Scalar,rows> >(u);
}

// -------------------------- multiplication and division by scalar -------------

template<typename Scalar, __size_type rows>
inline Vector<Scalar,rows> operator*(const Scalar&s, const ColumnVector<Scalar,rows>&u){
  return multiplyObjectScalar< Vector<Scalar,rows> > (u,s);
}

template<typename Scalar, __size_type rows>
inline Vector<Scalar,rows> operator*(const ColumnVector<Scalar,rows>&u, const Scalar&s){
  return multiplyObjectScalar< Vector<Scalar,rows> > (u,s);
}

template<typename Scalar, __size_type rows>
inline Vector<Scalar,rows> operator/(const ColumnVector<Scalar,rows>&u, const Scalar&s){
  return divideObjectScalar< Vector<Scalar,rows> > (u,s);
}

template<typename Scalar, __size_type rows>
inline ColumnVector<Scalar,rows>& ColumnVector<Scalar,rows>::operator*=(const Scalar& s){
  return multiplyAssignObjectScalar(*this,s);
}

template<typename Scalar, __size_type rows>
ColumnVector<Scalar,rows>& ColumnVector<Scalar,rows>::operator/=(const Scalar& s){
  return divideAssignObjectScalar(*this,s);
}

// -------------------------------------- assignments ---------------------------------------

template<typename Scalar, __size_type rows>
inline ColumnVector<Scalar,rows>& ColumnVector<Scalar,rows>::operator+=(const ColumnVector& v){
  return addAssignObjectObject(*this,v);
}

template<typename Scalar, __size_type rows>
inline ColumnVector<Scalar,rows>& ColumnVector<Scalar,rows>::operator+=(const Vector<Scalar,rows>& v){
  return addAssignObjectObject(*this,v);
}

template<typename Scalar, __size_type rows>
inline ColumnVector<Scalar,rows>& ColumnVector<Scalar,rows>::operator-=(const ColumnVector& v){
  return subtractAssignObjectObject(*this,v);
}

template<typename Scalar, __size_type rows>
inline ColumnVector<Scalar,rows>& ColumnVector<Scalar,rows>::operator-=(const Vector<Scalar,rows>& v){
  return subtractAssignObjectObject(*this,v);
}


template<typename Scalar, __size_type rows>
void ColumnVector<Scalar,rows>::clear(){
  capd::vectalg::clear(*this);
}



/// It serializes a matrix - gives text reprezentation which can be compiled
template<typename Scalar, __size_type rows>
std::string cppReprezentation(const ColumnVector<Scalar,rows> & A, const std::string& varName,
			      const std::string& typeName);

}} // namespace capd::vectalg

#endif // _CAPD_VECTALG_COLUMNVECTOR_H_

/// @}
