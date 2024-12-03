/// @addtogroup diffAlgebra
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file HomogenousPolynomial.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DIFFALGEBRA_HOMOGENOUS_POLYNOMIAL_H_
#define _CAPD_DIFFALGEBRA_HOMOGENOUS_POLYNOMIAL_H_

#include "capd/basicalg/factrial.h"
#include "capd/vectalg/Multiindex.h"

namespace capd{
namespace diffAlgebra{

/// Class HomogenousPolynomial provides indexing and some algorithms for multivariate homogenous polynomials.
/// It does not store the coefficients of the polynomial. It is assumed that the memory is already allocated in a continuous block.

template<class Scalar>
class HomogenousPolynomial
{
public:
  typedef Scalar ScalarType;
  typedef ScalarType* iterator;
  typedef ScalarType* const_iterator;;
  typedef ScalarType& reference;
  typedef const ScalarType& const_reference;
  typedef capd::vectalg::Multipointer Multipointer;
  typedef capd::vectalg::Multiindex Multiindex;

  HomogenousPolynomial(ScalarType* coefficients, int dim, int degree)
    : m_coefficients(coefficients), m_dim(dim), m_degree(degree)
  {}

  // general indexing
  reference operator()(const Multipointer&); //< returns reference to coefficient corresponding to a multipointer
  reference operator()(const Multiindex&);   //< returns reference to coefficient corresponding to a multiindex
  const_reference operator()(const Multipointer&) const; //< returns const reference to coefficient corresponding to a multipointer
  const_reference operator()(const Multiindex&) const; //< returns const reference to coefficient corresponding to a multiindex

  // fast indexing for C^1, C^2 and C^3 algorithms
  reference operator()(int i); //< returns reference to a coefficient in linear part, i.e. \f d/dx_i \f
  reference operator()(int i, int j); //< returns reference to a coefficient in second order part, i.e. \f d^2/dx_idx_j \f
  reference operator()(int i, int j, int c); //< returns reference to a coefficient in third order part, i.e. \f d^3/dx_idx_jdx_c \f
  const_reference operator()(int i) const; //< returns const reference to a coefficient in linear part, i.e. \f d/dx_i \f
  const_reference operator()(int i, int j) const; //< returns const reference to a coefficient in second order part, i.e. \f d^2/dx_idx_j \f
  const_reference operator()(int i, int j, int c) const; //< returns const reference to a coefficient in third order part, i.e. \f d^3/dx_idx_jdx_c \f

  iterator begin(){ return m_coefficients; } //< iterator selection. Returns iterator to the first element.
  iterator end(){ return m_coefficients + newton(m_dim+m_degree-1,m_degree); } //< iterator selection. Returns iterator to the element after last element.
  const_iterator begin() const { return m_coefficients; } //< const_iterator selection. Returns const_iterator to the first element.
  const_iterator end() const { return m_coefficients + newton(m_dim+m_degree-1,m_degree); } //< const_iterator selection. Returns const_iterator to the element after last element.
  reference operator[](int i){ return m_coefficients[i]; } //< access to an element by its direct position. Can be used for fast manipulation on all coefficients (like multiplication by constant).
  const_reference operator[](int i) const { return m_coefficients[i]; } //< const_reference to an element by its direct position. Can be used for fast manipulation on all coefficients (like multiplication by constant).

protected:
  ScalarType* m_coefficients; //< a pointer to externally allocated memory.
  int m_dim; //< number of variables in the polynomial.
  int m_degree; //< degree of homomgenous polynomial.
};

}} // namespace capd::diffAlgebra

#endif // _CAPD_DIFFALGEBRA_HOMOGENOUS_POLYNOMIAL_H_

/// @}
