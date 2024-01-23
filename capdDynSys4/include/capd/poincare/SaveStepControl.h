/////////////////////////////////////////////////////////////////////////////
/// @file BasicPoincareMap.h
///
/// @author Daniel Wilczak
/// @author Tomasz Kapela
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2008 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_SAVE_STEP_CONTROL_MAP_H_
#define _CAPD_SAVE_STEP_CONTROL_MAP_H_

namespace capd{
namespace poincare{
/// @addtogroup poincare 
/// @{

////////////////////////////////////////////////////////////////////////
///
/// it saves step and step control settings on construction
/// and restores them on destruction.
///
////////////////////////////////////////////////////////////////////////
template <typename DS>
class SaveStepControl{
public:
  typedef typename DS::ScalarType ScalarType;
  SaveStepControl(DS & ds)
  : m_ds(ds),
    //m_step(ds.getStep()),
    m_maxStep(ds.getMaxStep()),
    m_stepChangeAllowance(ds.isStepChangeAllowed()){
  }
  ~SaveStepControl(){
    //m_ds.setStep(m_step);
    m_ds.setMaxStep(m_maxStep);
    m_ds.onOffStepControl(m_stepChangeAllowance);
  }
protected:
  DS & m_ds;
  //ScalarType m_step;
  ScalarType m_maxStep;
  bool m_stepChangeAllowance;
};


/// @}
}} // namespace capd::poincare



#endif  // _CAPD_SAVE_STEP_CONTROL_MAP_H_

