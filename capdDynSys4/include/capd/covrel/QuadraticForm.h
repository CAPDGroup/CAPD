/////////////////////////////////////////////////////////////////////////////
/// @file QuadraticForm.h
///
/// @author Tomasz Kapela  @date 2009-08-03
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2016 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_COVREL_QUADRATICFORM_H_
#define _CAPD_COVREL_QUADRATICFORM_H_

namespace capd {
namespace covrel {
/// @addtogroup covrel
/// @{

template <typename MatrixT>
class QuadraticForm{
public:
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  /// General constructor
  QuadraticForm(const MatrixType & Q): m_Q(Q){
  }
  /// Constructor for dimension 2
  /// Q((x,y)) = ax^2 + by^2
  QuadraticForm(
      ScalarType a = capd::TypeTraits<ScalarType>::one(),
      ScalarType b = -capd::TypeTraits<ScalarType>::one()
  ) : m_Q(2,2) {
    m_Q[0][0] = a;              m_Q[0][1]=ScalarType(0.);
    m_Q[1][0] = ScalarType(0.); m_Q[1][1] = b;
  }
  /// Value of quadratic form
  ScalarType operator()(const VectorType & x){
	  return x*(m_Q*x);
  }
  // Derivative
  VectorType operator[](const VectorType & x){
    return Transpose(m_Q)*x+m_Q*x;
  }
  const MatrixType & getQ() const{
    return m_Q;
  }
  virtual std::string show(void) const{
    std::ostringstream str;
    str << "\n Quadratic form : " << m_Q;
    return str.str();
  }
protected:
  MatrixType m_Q;
};

/// @}
}} // namespace capd::covrel
#endif
