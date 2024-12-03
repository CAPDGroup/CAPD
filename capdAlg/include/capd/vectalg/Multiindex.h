/// @addtogroup vectalg
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file Multiindex.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_VECTALG_MULTIINDEX_H_
#define _CAPD_VECTALG_MULTIINDEX_H_

#include <stdexcept>
#include <vector>
#include "capd/basicalg/factrial.h"
#include "capd/basicalg/TypeTraits.h"
#include "capd/vectalg/Vector.h"

namespace capd{
namespace vectalg{


class Multiindex;


/// Multipointer always contains nondecreasing list of indexes of variables.
///
/// For partial derivatives they denote variables with respect to which we differentiate.
/// For example, a Multipointer mp=(0,0,2,3) corresponds to partial derivative
/// \f$ \frac{\partial^4}{\partial x_0^2 \partial x_2 \partial x_3} \f$
///
/// For power series multipointer denotes variables that form monomial.
/// e.g.Multipointer mp=(0,0,2,3) correspond to monomial \f$ x_0^2x_2x_3\f$
///
class Multipointer : public Vector<int,0>
{
public:
  typedef Vector<int,0>::iterator iterator;
  typedef Vector<int,0>::const_iterator const_iterator;
  typedef Vector<int,0>::size_type size_type;

  Multipointer(void){}
  explicit Multipointer(size_type _dimension) : Vector<int,0>(_dimension){}
  explicit Multipointer(const Multiindex&);
  Multipointer(const Multipointer& p) : Vector<int,0>(p) {}
  Multipointer(size_type dim, int data[]) : Vector<int,0>(dim,data){}
  Multipointer(size_type dim,bool) : Vector<int,0>(dim,true){}
  Multipointer& operator=(const Multipointer& v) {
    Vector<int,0>::operator= ( static_cast< const Vector<int,0> &>(v));
    return *this;
  }

  Multipointer(Multipointer&& v) : Vector<int,0>(std::move(v)) {}
  Multipointer & operator=(Multipointer && v) {
    Vector<int,0>::operator= ( static_cast< Vector<int,0> &&>(v));
    return *this;
  }
  Multipointer(std::initializer_list<int> l) : Vector<int,0>(l.size(), false) {
    std::copy(l.begin(), l.end(), this->begin());
  }

  /// order of derivative
  inline size_type module() const{
    return dimension();
  }
  bool hasNext(size_type dim);
  long factorial() const;
  Multipointer subMultipointer(const Multipointer& mp) const;

  typedef std::vector<Multipointer> MultipointersVector;
  typedef std::vector<MultipointersVector> IndicesSet;
  static const IndicesSet& generateList(size_type p, size_type k);

  size_type index(size_type dimension, size_type maxDegree) const;
  size_type index(size_type dimension, size_type maxDegree, const Multipointer& sub) const;

private:
  static std::vector<IndicesSet> indices;
  static void computeNextLevel();
  static size_type maxKnownLevel;
  inline static IndicesSet& getList(size_type n, size_type k)
  {
    return indices[n*(n-1)/2+k-1];
  }
};

// -------------------------------------------------------------------------------

/// For a Multiindex mi, mi[p] is a number of differentiation with respect to i-th variable.
/// For example, a Multipointer mp=(0,0,2,3) in 5-dimensional space corresponds to
/// the Multiindex mi=(2,0,1,1,0).
/// Hence, Multiindex agrees with standard notation and it contains an additional information
/// about the dimension of the domain of the function.
///
/// For polynomial:  Multiindex stores exponents of a given monomial.
/// e.g. monomial \f$ x^2 z^3 \f$ of 4 variables (x,y,z,t) has multiindex (2,0,3,0)
///
class Multiindex : public Vector<int,0>
{
public:
  typedef Vector<int,0>::iterator iterator;
  typedef Vector<int,0>::const_iterator const_iterator;
  typedef std::vector<Multiindex> MultiindexVector;
  typedef std::vector<MultiindexVector> IndicesSet;
  typedef Vector<int,0>::size_type size_type;

  Multiindex(void){}
  explicit Multiindex(size_type _dimension) : Vector<int,0>(_dimension){}
  Multiindex(size_type _dimension, const Multipointer&);
  Multiindex(size_type dim, int data[]) : Vector<int,0>(dim,data){}
  Multiindex(size_type dim,bool) : Vector<int,0>(dim,true){}
  Multiindex(const Multiindex& v) : Vector<int,0>(v) {}
  Multiindex& operator=(const Multiindex & v) {
    Vector<int,0>::operator= ( static_cast< const Vector<int,0> &>(v));
    return *this;
  }


  Multiindex(Multiindex&& v) : Vector<int,0>(std::move(v)) {}
  Multiindex & operator=(Multiindex && v) {
    Vector<int,0>::operator= ( static_cast< Vector<int,0> &&>(v));
    return *this;
  }
  Multiindex(std::initializer_list<int> l) : Vector<int,0>(l.size(), false) {
    std::copy(l.begin(), l.end(), this->begin());
  }

  bool hasNext();
  bool hasNext(int* a, int* b) const;
  bool hasNext(int* a, int* b, size_type j) const;

  size_type module() const;          ///< returns sum of multiindex coordinates (i.e. degree of monomial)
  long factorial() const;      ///< for multiindex (a,b,..,n) returns a!b!...n!
  static void generateList(size_type n, size_type k, IndicesSet& result);
  size_type index(size_type maxDegree) const; ///< computes index in the array that corresponds to this multiindex
};

// -------------------------------------------------------------------------------

Multipointer sumMultipointers(const Multipointer&, const Multipointer&);

/// returns new multipointer which is multiindex mp with index added in correct place
Multipointer addIndex(const Multipointer & mp, int index);

/// appends index to the end of multipointer mp
Multipointer push_back(const Multipointer & mp, int index);

/**
  checks if muiltipointer contains index
  @param[in] - mp multipointer that is checked
  @param[in] - index index to be found
  @return true if mp contains index
*/
bool hasIndex(const  Multipointer & mp, int index);

/**
  returns the number of occurences of index in the multipointer
  @param[in] mp - multipointer that is inspected
  @param[in] index - the counted index
  @return number of occureneces o of index in mp
*/
int indexCount(const Multipointer &mp, int index);

/**
  returns a multipointer with removed index
  @param[in] mp - multipointer from which an index is removed
  @param[in] index -index to be removed
  @returns a multipointer with removed index
 */
Multipointer removeIndex(const Multipointer & mp, int index);

}} // namespace capd::vectalg

// -------------------------------------------------------------------------------

///
/// It computes v^m where v is a vector and m is a multiindex
///

template<typename VectorType>
typename VectorType::ScalarType power(const VectorType& v, const capd::vectalg::Multiindex& m)
{
  using namespace capd::vectalg;
  if(v.dimension()!=m.dimension())
    throw std::runtime_error("power(vector,multiindex) error: different dimensions of vector and multiindex");
  typename VectorType::ScalarType result=1;
  typename VectorType::const_iterator b=v.begin(), e=v.end();
  typename Multiindex::const_iterator i=m.begin();
  while(b!=e)
  {
    result *= power(*b,*i);
    ++b;
    ++i;
  }
  return result;
}

// -------------------------------------------------------------------------------

template<typename VectorType>
typename VectorType::ScalarType power(const VectorType& v, const capd::vectalg::Multipointer& m)
{
  using namespace capd::vectalg;
  typedef typename VectorType::ScalarType ScalarType;
  ScalarType result = capd::TypeTraits<ScalarType>::one();
  typename Multipointer::const_iterator b=m.begin(), e=m.end();
  while(b!=e)
  {
    typename Multipointer::const_iterator temp=b;
    int p = *b;
    do{
      ++b;
    }while(b!=e && *b==p);
    size_t n = b-temp;
    result *= power(v[p],(int)n);
  }
  return result;
}

#endif // _CAPD_VECTALG_MULTIINDEX_H_

/// @}
