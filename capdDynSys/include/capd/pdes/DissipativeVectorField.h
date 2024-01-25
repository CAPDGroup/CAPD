/// @addtogroup pdes
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file DissipativeVectorField.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2008-2016 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_PDES_DISSIPATIVEVECTORFIELD_H_
#define _CAPD_PDES_DISSIPATIVEVECTORFIELD_H_

#include <vector>
#include "capd/intervals/lib.h"
#include "capd/vectalg/lib.h"

namespace capd {
namespace pdes {

/**
 * The class provides a common interface for a dissipative vector required by the class PdeSolver.
 * Each particular system should inherit and implement abstract methods.
*/
template<class SeriesT>
class DissipativeVectorField {
public:
  // typedef SeriesT SeriesType;
  typedef capd::interval ScalarType;
  typedef SeriesT VectorType;
  typedef capd::IMatrix MatrixType;
  typedef MatrixType::RowVectorType FiniteVectorType;

  typedef capd::IVector::size_type size_type;
  typedef std::vector<VectorType> VectorArray;
  typedef std::vector<VectorArray> MatrixArray;

  virtual VectorType operator()(ScalarType h, const VectorType& v) = 0;
  virtual VectorType operator()(ScalarType h, const VectorType& v, MatrixType& DF) = 0;
  virtual MatrixType derivative(ScalarType h, const VectorType& v) = 0;

  /// This function should compute ODE coefficients up to given order at the set a
  virtual void computeODECoefficients(VectorArray& a, size_type order) = 0;

  /// Given coefficients for C^0 part, it computes Taylor coefficient for one column of variational equation
  /// Initial condition for this column is 'c'.
  virtual void computeODECoefficients(const VectorArray& a, VectorArray& c, size_type order) = 0;

  /// This function should compute ODE coefficients up to given order at the set a
  /// Moreover, block derivative of first group of variables has to be computed.
  void computeODECoefficients(VectorArray& a, MatrixArray& J, size_type p, size_type numberOfColumns){
    this->computeODECoefficients(a,p);
    for(size_type c=0;c<numberOfColumns;++c)
      this->computeODECoefficients(a,J[c],p);
  }

  /// This function should refine the tail so that the vector field is pointing inwards the tail.
  virtual void makeSelfConsistentBound(VectorArray& a) = 0;

  /// This function should refine the tail so that the vector field for variational equation is pointing inwards the tail.
  /// Here we assume two different initial conditions J1=(Id,0) and J2=(0,something)
  virtual void makeSelfConsistentBound(VectorArray& a, MatrixArray& J1, MatrixArray& J2, size_type numberOfColumns) = 0;

  virtual size_type dimension() const = 0;
  virtual size_type firstDissipativeIndex() const = 0;

  // these functions update tail of a vector after time step h using differential inequality and under assumptions that isolation is satisfied
  /// Update tail for C^0 part using linear differential inequality
  virtual void updateTail(VectorType& x, const VectorArray& enc, ScalarType h) const = 0;
  /// Update tail for two C^1 blocks using linear differential inequality
  virtual void updateTail(VectorArray& DyxId, VectorArray& Dyx, const MatrixArray& Enc, const MatrixArray& DyxEnc, ScalarType h) const = 0;

  /// This function should compute a matrix M such that
  /// M_ii is logarithmic norm of the diagonal block
  /// M_ij is a norm of ij block
  /// The infinite dimensional space is split onto m+1 blocks
  virtual MatrixType blockNorms(const VectorType& a, size_type m) const = 0;
};

}} // namespace capd::pdes


#endif // _CAPD_PDES_DISSIPATIVEVECTORFIELD_H_

/// @}
