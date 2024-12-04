/////////////////////////////////////////////////////////////////////////////
/// @file TimeRange.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DIFFALGEBRA_TIMERANGE_H_
#define _CAPD_DIFFALGEBRA_TIMERANGE_H_

namespace capd{
namespace diffAlgebra{
/// @addtogroup diffAlgebra
/// @{

///////////////////////////////////////////////////////////////////////////////////
/// TimeRange is a base class for all types of sets. It stores the current time in ODE case
/// or number of iterations in the case of discrete dynamical system case.
///
//////////////////////////////////////////////////////////////////////////////////

template<class IntervalT>
class TimeRange{
public:
  typedef IntervalT ScalarType;

  explicit TimeRange(const ScalarType& t)
    : m_currentTime(t)
  {}

  const ScalarType getCurrentTime() const{
    return m_currentTime;
  }
  ScalarType& refCurrentTime() {
    return m_currentTime;
  }
  void setCurrentTime(const ScalarType& t){
    m_currentTime = t;
  }
protected:
  ScalarType m_currentTime;
};

/// @}
}} // the end of the namespace capd::diffAlgebra

#endif // _CAPD_DIFFALGEBRA_TIMERANGE_H_

