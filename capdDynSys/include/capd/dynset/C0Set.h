/////////////////////////////////////////////////////////////////////////////
/// @file C0Set.h
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_C0SET_H_
#define _CAPD_DYNSET_C0SET_H_

#include <string>
#include <sstream>
#include <stdexcept>
#include "capd/basicalg/factrial.h"
#include "capd/diffAlgebra/TimeRange.h"
#include "capd/dynset/EnclosureHolder.h"
#include "capd/dynsys/DynSys.h"
#include "capd/dynset/SetTraits.h"

namespace capd{
/**
 * Various set representations that can be moved with dynamical systems
 */
namespace dynset{
/// @addtogroup dynset
/// @{

/**
 * Common interface of all sets that stores only C0 information (set position).
 */
template<typename MatrixT>
class C0Set : public capd::diffAlgebra::TimeRange<typename MatrixT::ScalarType>,
              public capd::dynset::C0EnclosureHolder<typename MatrixT::RowVectorType>
{
public:
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;
  typedef C0Set SetType;                                   ///<  defines class of given set (C0, C1, C2, ...)
  typedef capd::dynsys::DynSys<MatrixType> DynSysType;     ///< defines minimal regularity of the dynamical system

  /// constructor, initializes default enclosure and initial time
  C0Set(const VectorType& x, const VectorType& enc, const ScalarType& t)
      : capd::diffAlgebra::TimeRange<ScalarType>(t),
        capd::dynset::C0EnclosureHolder<VectorType>(x,enc)
  {}

  /// destructor
  virtual ~C0Set (void) {}
    /// computes image of the set after one step/iterate of the dynamical system
  virtual void move(DynSysType& dynsys) = 0;
  const static size_type degree() { return 0; }
};

/**
 * Specialization of Traits class
 */

template<class MatrixT>
struct SetTraits< C0Set<MatrixT> >{
	const static bool isC0Set=true;
	const static bool isC1Set=false;
	const static bool isC2Set=false;
	const static bool isCnSet=false;
};

/**
 * This is an algorithm for computation of polynomial part of Hermite-Obreshkov interpolation
 * Used in all HO-sets
 *
 * Assumptions:
 * @param ds dynamical system. We assume that it has computed coefficients for the main part and variational equations
 * @param p,q two numners such that p+q=order
 * @param[out] psi computed polynomial interpolation
 * @param[out] Dpsi derivative of psi wrt
 *
 * NOTE: In the present implementation order cannot be bigger than 33. This is due to the capacity of long type when compute binomial.
 */

template<class DS, class V, class M>
void computePsi(DS& ds, int p, int q, const typename DS::ScalarType& h, V& psi, M& Dpsi)
{
  int dimension = psi.dimension();

  int i,j,k;
  unsigned long denominator = binomial(p+q,p);
  for(i=0;i<dimension;++i)
  {
    psi[i] = ds.centerCoefficient(i,p)/denominator;
    for(j=0;j<dimension;++j)
      Dpsi(i+1,j+1) = ds.coefficient(i,j,p)/denominator;
  }

  for(k=p-1;k>=0;k--)
  {
    for(i=0;i<dimension;++i)
    {
      psi[i] = psi[i]*h + (unsigned long)binomial(p+q-k,q)*ds.centerCoefficient(i,k)/denominator;
      for(j=0;j<dimension;++j)
        Dpsi(i+1,j+1) = Dpsi(i+1,j+1)*h + (unsigned long)binomial(p+q-k,q)*ds.coefficient(i,j,k)/denominator;
    }
  }
} // end of computePsi

/// @}
}} //namespace capd::dynset

#endif // _CAPD_DYNSET_C0SET_H_
