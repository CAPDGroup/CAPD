/// @addtogroup vectalg
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file Container.h
///
/// This file provides a template class Container together with suitable iterators
/// The container has fixed size specified by a template argument 'capacity'
///
/// Also a specialization of this class for capacity=0 is defined
/// In that case objects in this container are allocated on storage instead of stack
///
/// This class is used as a container for vectors, matrices and higher order containers
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.


#ifndef _CAPD_VECTALG_CONTAINER_H_
#define _CAPD_VECTALG_CONTAINER_H_

#include <stdexcept>
#include <cstdlib> // for size_t
#include "capd/settings/compilerSetting.h"
#include "capd/auxil/Dll.h"

namespace capd{
namespace vectalg{

typedef unsigned __size_type;
typedef int __difference_type; // must be signed integral type

/// class Container together with suitable iterators
/// The container has fixed size specified by a template argument 'capacity'
///
/// This class is used as a container for vectors, matrices and higher order containers
template<typename Scalar, __size_type capacity>
class Container
{
public:
  typedef Scalar ScalarType;
  typedef __size_type size_type;
  typedef __difference_type difference_type;
  typedef ScalarType* iterator;
  typedef const ScalarType* const_iterator;
  typedef std::reverse_iterator<iterator> reverse_iterator;
  typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

  Container();
  explicit Container(size_type);
  Container(size_type,bool); // it does not insert zeros

  iterator begin();
  iterator end();
  const_iterator begin() const;
  const_iterator end() const;

  reverse_iterator rbegin();
  reverse_iterator rend();
  const_reverse_iterator rbegin() const;
  const_reverse_iterator rend() const;

  /* Container& operator=(const Container&);*/
  void resize(size_type newCapacity);

  ScalarType& operator[](size_type);
  const ScalarType& operator[](size_type) const;
  ScalarType& operator()(size_type);
  const ScalarType& operator()(size_type) const;
  friend void swap(Container<Scalar,capacity>& A_c1, Container<Scalar,capacity>& A_c2) noexcept
  {
    iterator b=A_c1.begin(), e=A_c1.end();
    iterator i = A_c2.begin();
    while(b!=e)
    {
      std::swap(*b,*i);
      ++b;
      ++i;
    }
  }
  void clear();
// memory allocation
  static size_type size() {return capacity;}

protected:
  ScalarType data[capacity];
};

/// Specialization for capacity=0
/// This container allocates objects on a storage
template<typename Scalar>
class Container<Scalar,0>
{
public:
  typedef Scalar ScalarType;
  typedef __size_type size_type;
  typedef __difference_type difference_type;
  typedef ScalarType* iterator;
  typedef const ScalarType* const_iterator;
  typedef std::reverse_iterator<iterator> reverse_iterator;
  typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

  Container& operator=(Container&&) noexcept;
  Container(Container&&) noexcept;

  Container();
  explicit Container(size_type);
  Container(size_type,bool); // it does not insert zeros
  Container(const Container&);
  ~Container();

  iterator begin();
  iterator end();
  const_iterator begin() const;
  const_iterator end() const;

  reverse_iterator rbegin();
  reverse_iterator rend();
  const_reverse_iterator rbegin() const;
  const_reverse_iterator rend() const;

  Container& operator=(const Container&);
  void resize(size_type);

  ScalarType& operator[](size_type);
  const ScalarType& operator[](size_type) const;
  ScalarType& operator()(size_type);
  const ScalarType& operator()(size_type) const;
  friend void swap(Container<Scalar,0>& A_c1, Container<Scalar,0>& A_c2) noexcept
  {
     std::swap(A_c1.data,A_c2.data);
     std::swap(A_c1.capacity,A_c2.capacity);
  }

  size_type size() const { return capacity; }
  bool empty() const { return capacity == 0; }

  void clear();

protected:
  ScalarType *data = nullptr;
  size_type capacity=0u;
};

// ---- inline definitions for Containers ----------------- //

template<typename Scalar>
inline Container<Scalar,0>::Container() : data(0), capacity(0u)
{
}
// --------------- iterator selection --------------------- //

template<typename Scalar, __size_type capacity>
inline typename Container<Scalar,capacity>::iterator Container<Scalar,capacity>::begin()
{
  return iterator(data);
}

template<typename Scalar, __size_type capacity>
inline typename Container<Scalar,capacity>::iterator Container<Scalar,capacity>::end()
{
  return iterator(data+capacity);
}

template<typename Scalar, __size_type capacity>
inline typename Container<Scalar,capacity>::const_iterator Container<Scalar,capacity>::begin() const
{
  return const_iterator(data);
}

template<typename Scalar, __size_type capacity>
inline typename Container<Scalar,capacity>::const_iterator Container<Scalar,capacity>::end() const
{
  return const_iterator(data+capacity);
}

template<typename Scalar, __size_type capacity>
inline typename Container<Scalar,capacity>::reverse_iterator Container<Scalar,capacity>::rbegin()
{
  return reverse_iterator(end());
}

template<typename Scalar, __size_type capacity>
inline typename Container<Scalar,capacity>::reverse_iterator Container<Scalar,capacity>::rend()
{
  return reverse_iterator(begin());
}

template<typename Scalar, __size_type capacity>
inline typename Container<Scalar,capacity>::const_reverse_iterator Container<Scalar,capacity>::rbegin() const
{
  return const_reverse_iterator(end());
}

template<typename Scalar, __size_type capacity>
inline typename Container<Scalar,capacity>::const_reverse_iterator Container<Scalar,capacity>::rend() const
{
  return const_reverse_iterator(begin());
}

template<typename Scalar>
inline typename Container<Scalar,0>::iterator Container<Scalar,0>::begin()
{
  return iterator(data);
}

template<typename Scalar>
inline typename Container<Scalar,0>::iterator Container<Scalar,0>::end()
{
  return iterator(data+capacity);
}

template<typename Scalar>
inline typename Container<Scalar,0>::const_iterator Container<Scalar,0>::begin() const
{
  return const_iterator(data);
}

template<typename Scalar>
inline typename Container<Scalar,0>::const_iterator Container<Scalar,0>::end() const
{
  return const_iterator(data+capacity);
}

template<typename Scalar>
inline typename Container<Scalar, 0>::reverse_iterator Container<Scalar, 0>::rbegin()
{
  return reverse_iterator(end());
}

template<typename Scalar>
inline typename Container<Scalar, 0>::reverse_iterator Container<Scalar, 0>::rend()
{
  return reverse_iterator(begin());
}

template<typename Scalar>
inline typename Container<Scalar, 0>::const_reverse_iterator Container<Scalar, 0>::rbegin() const
{
  return const_reverse_iterator(end());
}

template<typename Scalar>
inline typename Container<Scalar, 0>::const_reverse_iterator Container<Scalar, 0>::rend() const
{
  return const_reverse_iterator(begin());
}

// ------------------------- indexing ------------------------ //

template<typename Scalar, __size_type capacity>
inline Scalar& Container<Scalar,capacity>::operator[] (size_type i)
{
  return data[i];
}

template<typename Scalar, __size_type capacity>
inline const Scalar& Container<Scalar,capacity>::operator[] (size_type i) const
{
  return data[i];
}

template<typename Scalar, __size_type capacity>
inline Scalar& Container<Scalar,capacity>::operator() (size_type i)
{
  return data[i-1];
}

template<typename Scalar, __size_type capacity>
inline const Scalar& Container<Scalar,capacity>::operator() (size_type i) const
{
  return data[i-1];
}

template<typename Scalar>
inline Scalar& Container<Scalar,0>::operator[] (size_type i)
{
  return data[i];
}

template<typename Scalar>
inline const Scalar& Container<Scalar,0>::operator[] (size_type i) const
{
  return data[i];
}

template<typename Scalar>
inline Scalar& Container<Scalar,0>::operator() (size_type i)
{
  return data[i-1];
}

template<typename Scalar>
inline const Scalar& Container<Scalar,0>::operator() (size_type i) const
{
  return data[i-1];
}

// ------------ constructor - desctructor --------------------

template<typename Scalar, __size_type capacity>
inline Container<Scalar,capacity>::Container(size_type,bool)
{}

template<typename Scalar>
inline Container<Scalar,0>::Container(size_type a_capacity,bool) : capacity(a_capacity)
{
  data = new ScalarType[capacity];
}

template<typename Scalar>
inline Container<Scalar,0>::~Container()
{
  if(data) delete [] data;
}

}} // namespace capd::vectalg

#endif // _CAPD_VECTALG_CONTAINER_H_

/// @}
