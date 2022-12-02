/////////////////////////////////////////////////////////////////////////////
/// @file Vector.h
///
/// This file provides a template class Vector together with typical
/// algebraic operations. Most of them are defined as generic algorithms
/// in files 'commonOperations.h' and 'commonOperations.hpp'
/// For inline definitions of operators see to file
/// Vector_inline.h included at the end of this file.
///
/// The class uses class 'Container' as a container for coefficients
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.


#ifndef _CAPD_VECTALG_VECTOR_H_
#define _CAPD_VECTALG_VECTOR_H_

#include <iostream>
#include <cstdlib>
#include <vector>
#include "capd/basicalg/minmax.h"
#include "capd/basicalg/power.h"
#include "capd/vectalg/Container.h"
#include "capd/settings/compilerSetting.h"

namespace capd{
namespace vectalg{

/// @addtogroup vectalg
/// @{


template<typename Scalar, __size_type dim>
class Vector;

template<typename Scalar, __size_type dim>
std::ostream& operator<<(std::ostream& out, const Vector<Scalar,dim>& v);

template<typename Scalar, __size_type dim>
std::istream& operator>>(std::istream& inp, Vector<Scalar,dim>& v);


//########################### Vector template ################################//

template<typename Scalar, __size_type dim>
class Vector : public Container<Scalar,dim>
{
public:
  typedef Scalar ScalarType;
  typedef Container<Scalar,dim> ContainerType;
  typedef typename ContainerType::iterator iterator;
  typedef typename ContainerType::const_iterator const_iterator;
  typedef Vector<Scalar,dim> VectorType;
  typedef __size_type size_type;
  typedef __difference_type difference_type;

  Vector(void){}
  explicit Vector(size_type a_dim); // for compatibility with heap-allocated specialization
  Vector(const Scalar& x, const Scalar& y,const Scalar& z);  //obsolete, requires dimension=3
  Vector(size_type,const ScalarType[]);
  explicit Vector(const char data[]);
  explicit Vector(const std::string & data);
  Vector(const Vector&);
  template<typename S,
          typename std::enable_if<std::is_convertible<S, ScalarType>::value && !std::is_same<S, ScalarType>::value, int>::type = 0 >
  explicit Vector(const Vector<S,dim>&);

  Vector(size_type,bool); // it does not insert zeros

  template<__size_type dataDim>
  Vector(const Scalar (&data)[dataDim]);

  template<typename Iterator>
  Vector(Iterator begin, Iterator end);

  Vector(Vector&& v) = default;
  Vector & operator=(Vector && v)= default;
  Vector(std::initializer_list<ScalarType> l);

  // assignments - vectors
  Vector& operator=  (const Vector& v);      //<assign a vector
  Vector& operator+= (const Vector& v);      //<add and assign a vector
  Vector& operator-= (const Vector& v);      //<subtract and assign a vector

  // assignments - Scalars
  Vector& operator=  (const Scalar& s);       //<assign a Scalar to each coordinate
  Vector& operator+= (const Scalar& s);       //<component-wise increase by a Scalar
  Vector& operator-= (const Scalar& s);       //<component-wise decrease by a Scalar
  Vector& operator*= (const Scalar& s);       //<scale by multiplying by Scalar
  Vector& operator/= (const Scalar& s);       //<scale by dividing by Scalar

  template<typename U>
  struct rebind {
      typedef Vector<U,dim> other;
  };

  size_type dimension() const {return ContainerType::size();}
  using ContainerType::begin;
  using ContainerType::end;
  using ContainerType::rbegin;
  using ContainerType::rend;
  using ContainerType::resize;
  using ContainerType::clear;

  // Euclidean norm
  ScalarType euclNorm(void) const;
  //if possible vector is normalized and true is returned. Otherwise false is returned.
  bool normalize(void);
  void sorting_permutation(typename rebind<int>::other& perm);

  const static size_type csDim = dim;
  static size_type degree() {return 0;} // required interface for DynSys
  static Vector* makeArray(size_type N, size_type _dim);
protected:
  using ContainerType::size;

}; // end of class Vector template

template<typename Vector>
std::string vectorToString( const Vector & v, int firstIndex = 0, int lastIndex = -1, int precision = -1);

template<typename Vector>
std::ostream & printVector(std::ostream & str, const Vector & v, int firstIndex = 0, int lastIndex = -1);

template<typename Scalar, __size_type dim>
inline std::ostream & print(std::ostream & str, const Vector<Scalar, dim> & v, int firstIndex = 0, int lastIndex = -1){
  return printVector(str, v, firstIndex, lastIndex);
}

/// It serializes a matrix - gives text reprezentation which can be compiled
template<typename Scalar, __size_type dim>
std::string cppReprezentation(const Vector<Scalar,dim> & A, const std::string& varName,
			      const std::string& typeName);

/// @}


}} // namespace capd::vectalg

#include "capd/vectalg/Vector_inline.h"

#endif // _CAPD_VECTALG_VECTOR_H_
