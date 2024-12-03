
/////////////////////////////////////////////////////////////////////////////
/// @file Jet.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include "capd/basicalg/factrial.h"
#include <vector>
#include "capd/diffAlgebra/CnContainer.h"
#include "capd/diffAlgebra/Hessian.h"

#ifndef _CAPD_DIFFALGEBRA_JET_H_
#define _CAPD_DIFFALGEBRA_JET_H_

namespace capd{
namespace diffAlgebra{
/// @addtogroup diffAlgebra
/// @{

/**
 * The class is used to store coefficients of a truncated power series to degree D
 * \f$ f:R^N->R^M \f$
 * Coefficients area assumed to be of a numeric type
 *
 * The base CnContainer is used to manage memory allocations, copying, etc.
*/

template<typename MatrixT, __size_type DEGREE>
class Jet : public CnContainer<
                  typename MatrixT::ScalarType,
                  MatrixT::ColumnVectorType::csDim,
                  MatrixT::RowVectorType::csDim,
                  DEGREE
                >
{
public:
  typedef MatrixT MatrixType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MatrixType::RefColumnVectorType RefVectorType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ColumnVectorType ImageVectorType;
  typedef capd::diffAlgebra::Hessian<ScalarType,ImageVectorType::csDim,VectorType::csDim> HessianType;

  typedef ScalarType* iterator;
  typedef const ScalarType* const_iterator;
  typedef CnContainer<ScalarType,ImageVectorType::csDim,VectorType::csDim,DEGREE> CnContainerType;
  typedef typename CnContainerType::BaseContainer ContainerType;
  typedef typename CnContainerType::Multipointer Multipointer;
  typedef typename CnContainerType::Multiindex Multiindex;
  typedef __size_type size_type;
  typedef __difference_type difference_type;

  Jet();
  Jet(size_type n, size_type d);
  Jet(size_type m, size_type n, size_type d);
  Jet(const Jet& v) = default;
  Jet(Jet&& v) : CnContainerType(static_cast< CnContainerType &&>(v)) {}
  Jet & operator=(const Jet & ) = default;
  Jet & operator=(Jet && v) {
    CnContainerType::operator= ( static_cast< CnContainerType &&>(v));
    return *this;
  }
  /// this operator returns a vector of partial derivatives, i.e. result[i] = d^{mp}f_i
  RefVectorType operator()(const Multipointer& mp) const;

  /// this operator returns a vector of partial derivatives, i.e. result[i] = d^{mp}f_i
  RefVectorType operator()(const Multiindex& mp) const;

  /// this operator returns a value of function, i.e. 0-order derivatives
  operator ImageVectorType() const;

  /// this operator returns value of function as a reference
  RefVectorType operator()(void) const;

  /// this operator returns first order derivatives as a matrix
  operator MatrixType() const;
  void setMatrix(const MatrixType&);

  /// this operator returns second order derivatives as a hessian object
  operator HessianType() const;

  // since CnCoeff represents a power series
// it may be used as a function

  ImageVectorType operator()(const VectorType&) const;
  MatrixType derivative(const VectorType& v) const;
  ImageVectorType operator()(const VectorType& v, MatrixType& matrix) const;

  template<typename U>
  struct rebind {
      typedef Jet<U,DEGREE> other;
  };

  /// returns string containing derivatives information
  std::string toString(int minFun = 0, int maxFun = -1,
                       int firstVariable = 0, int minDegree = 0, int maxDegree= -1,
                       int precision = -1) const;
  /// saves data to stream out with given precision
  /// (the default precision allows to load data without lose of precision)
  std::ostream & save(
      std::ostream & out,
      std::streamsize prec = capd::TypeTraits<ScalarType>::numberOfDigits() + 1
  ) const;
  friend std::ostream & operator << (std::ostream & out, const Jet & jet){
      return jet.save(out, 0);
  }

  /// loads cn data from given stream
  /// (expected format is exactly this provided by save)
  std::istream & load(std::istream & in);
  friend std::istream & operator >> (std::istream & in, Jet & jet){
    return jet.load(in);
  }

  using CnContainerType::operator();
  using CnContainerType::begin;
  using CnContainerType::end;
  using CnContainerType::dimension;
  using CnContainerType::degree;
  using CnContainerType::resize;

//protected:
  using CnContainerType::m_N;
  using CnContainerType::m_M;
}; // the end of class Jet


template<typename MatrixT, __size_type DEGREE>
inline bool operator== (const Jet<MatrixT,DEGREE>& x, const Jet<MatrixT,DEGREE>& y){
  return capd::vectalg::equal(x,y);
}

template<typename MatrixT, __size_type DEGREE>
inline bool operator!= (const Jet<MatrixT,DEGREE>& x, const Jet<MatrixT,DEGREE>& y){
  return capd::vectalg::notEqual(x,y);
}

/**
 * prints derivatives in the human readable form
 *
 * It prints all derivatives (df_i/dx_j^a)
 * \f$ \frac{ \partial^|a| f_i}{\partial x_j^a} \f$
 * where
 *   |a| = minRank, ..., maxRank
 *   i   = firstFun, ..., lastFun
 *   j   = firstVariable, ..., dimension.
 *
 *   a is a multiindex.
 */
template<typename MatrixT, __size_type DEGREE>
std::ostream & print(
    std::ostream & str,                           ///< output stream
    const Jet<MatrixT,DEGREE> & coeff,            ///< Cn coefficients
    int minDegree = 0,
    int maxDegree = -1,                           ///< default value is coeff.degree()
    int firstFun = 0,
    int lastFun = -1,                             ///< default value is coeff.imageDimension()
    int firstVariable = 0
);

// --------------------- inline definitions ------------------------

template<typename MatrixT, __size_type DEGREE>
inline Jet<MatrixT,DEGREE>::Jet()
    : CnContainerType(1,1,1)
{
  this->clear();
}

// -----------------------------------------------------------------

template<typename MatrixT, __size_type DEGREE>
inline Jet<MatrixT,DEGREE>::Jet(size_type dim, size_type degree)
  : CnContainerType(dim, dim, degree)
{
  this->clear();
}

// -----------------------------------------------------------------

template<typename MatrixT, __size_type DEGREE>
inline Jet<MatrixT,DEGREE>::Jet(size_type m, size_type n, size_type degree)
  : CnContainerType(m, n, degree)
{
  this->clear();
}

// -----------------------------------------------------------------

template<typename MatrixT, __size_type DEGREE>
inline typename Jet<MatrixT,DEGREE>::RefVectorType
Jet<MatrixT,DEGREE>::operator()(const Multipointer& mp) const
{
  return RefVectorType(
        begin(0,mp.dimension()) + mp.index(this->dimension(),this->degree()),
        binomial(this->dimension()+this->degree(),this->degree()),
        this->imageDimension()
      );
}

// -----------------------------------------------------------------

template<typename MatrixT, __size_type DEGREE>
inline typename Jet<MatrixT,DEGREE>::RefVectorType
Jet<MatrixT,DEGREE>::operator()(const Multiindex& mp) const
{
  return RefVectorType(
        begin(0,mp.module()) + mp.index(this->m_D),
        binomial(this->dimension()+this->degree(),this->degree()),
        this->imageDimension()
      );
}

// -----------------------------------------------------------------

template<typename MatrixT, __size_type DEGREE>
inline typename Jet<MatrixT,DEGREE>::RefVectorType
Jet<MatrixT,DEGREE>::operator()(void) const
{
  return RefVectorType(begin(),binomial(this->dimension()+this->degree(),this->dimension()),this->imageDimension());
}
/// @}
}} // namespace capd::diffAlgebra


// Since CnCoeff represents a power series (in fact all the partial derivatives of some function),
// there are obvious operators defined for the power series
template<typename MatrixT, capd::vectalg::__size_type DEGREE>
capd::diffAlgebra::Jet<MatrixT,DEGREE>
operator+(const capd::diffAlgebra::Jet<MatrixT,DEGREE>&, const capd::diffAlgebra::Jet<MatrixT,DEGREE>&);

template<typename MatrixT, capd::vectalg::__size_type DEGREE>
capd::diffAlgebra::Jet<MatrixT,DEGREE>
operator-(const capd::diffAlgebra::Jet<MatrixT,DEGREE>&, const capd::diffAlgebra::Jet<MatrixT,DEGREE>&);

template<typename MatrixT, capd::vectalg::__size_type DEGREE>
capd::diffAlgebra::Jet<MatrixT,DEGREE>
operator*(const MatrixT&, const capd::diffAlgebra::Jet<MatrixT,DEGREE>&);

template<typename MatrixT, capd::vectalg::__size_type DEGREE>
void substitutionPowerSeries(
      const capd::diffAlgebra::Jet<MatrixT,DEGREE>&,
      const capd::diffAlgebra::Jet<MatrixT,DEGREE>&,
      capd::diffAlgebra::Jet<MatrixT,DEGREE>& result,
      bool nonlinearOnly
   );

/**
 * This function computes inverse of power series that is close to identity, i.e. linear part is the Identity matrix.
 * @param[in] c power series that is inverted
 * @returns truncated power series of the inverse of c
 * @note linear part of c must be Identity matrix.
*/
template<class Jet>
Jet inverseSeriesCloseToIdentity(const Jet& c);

/**
 * This function computes inverse of a general power series.
 * @param[in] c power series that is inverted
 * @param[in] J inverse of linear part of series c.
 * @returns truncated power series of the inverse of c
 * @note In this procedure we do not specify how to compute inverse of linear part.
*/
template<class Jet>
Jet inversePowerSeries(const Jet& c, const typename Jet::MatrixType& J);

#endif // _CAPD_DIFFALGEBRA_JET_H_


