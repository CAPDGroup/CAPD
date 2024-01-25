/// @addtogroup diffAlgebra
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file CnTimeJet.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include <vector>
#include "capd/diffAlgebra/TimeRange.h"
#include "capd/diffAlgebra/Jet.h"
#include "capd/diffAlgebra/CoeffTraits.h"

#ifndef _CAPD_DIFFALGEBRA_CNTIMEJET_H_
#define _CAPD_DIFFALGEBRA_CNTIMEJET_H_

namespace capd{
namespace diffAlgebra{

/**
 * The class represent a jet of solution to a nonautomnomous ODE.
 * This is a truncated power series of degree D of a function
 * \f$ t \to R^n \f$
 * at given time
 *
 * The base class Jet is used to store coefficients of the jet.
 * The base class TimeRange is used to store current time.
*/

template<typename MatrixT, __size_type DEGREE>
class CnTimeJet : public capd::diffAlgebra::TimeRange<typename MatrixT::ScalarType>
{
public:
  typedef CnTimeJet SetType;
  typedef MatrixT MatrixType;
  typedef Jet<MatrixT,DEGREE> JetType;

  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MatrixType::RefColumnVectorType RefVectorType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ColumnVectorType ImageVectorType;
  typedef typename VectorType::size_type size_type;

  typedef typename JetType::iterator iterator;
  typedef typename JetType::const_iterator const_iterator;
  typedef typename JetType::Multipointer Multipointer;
  typedef typename JetType::Multiindex Multiindex;

  explicit CnTimeJet(JetType* jet)
    : capd::diffAlgebra::TimeRange<ScalarType>(TypeTraits<ScalarType>::zero()),
      m_jet(jet)
  {}

  // this operator returns a value of function, i.e. 0-order derivatives
  operator VectorType() const{
    return VectorType(*m_jet);
  }
  RefVectorType vector(){
    return (*m_jet)();
  }
  const RefVectorType vector() const{
    return (*m_jet)();
  }

  template<class DynSysT>
  void move(DynSysT& dynsys){
    dynsys(this->m_currentTime,*m_jet);
  }

  size_type degree() const{
    return m_jet->degree();
  }
protected:
  JetType* m_jet;
}; // the end of class CnCoeff

// --------------------- inline definitions ------------------------

template<typename MatrixT, __size_type DEGREE>
struct CoeffTraits<CnTimeJet<MatrixT,DEGREE> >{
	const static bool isC0Jet=false;
	const static bool isC1Jet=false;
	const static bool isC2Jet=false;
	const static bool isCnJet=true;
};
}} // namespace capd::diffAlgebra

#endif // _CAPD_DIFFALGEBRA_CNCOEFF_H_

/// @}
