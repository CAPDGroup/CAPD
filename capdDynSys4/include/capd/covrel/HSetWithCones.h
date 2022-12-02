//////////////////////////////////////////////////////////////////////////////
///
///  @file HSetWithCones.h
///  
///  @author kapela  @date   Aug 10, 2011
//////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2016 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_COVREL_HSETWITHCONES_H_
#define _CAPD_COVREL_HSETWITHCONES_H_
#include "capd/covrel/QuadraticForm.h"

namespace capd {
namespace covrel {
/// @addtogroup covrel
/// @{
  
template <typename HSetT, typename QFormT = QuadraticForm<typename HSetT::IMatrixType> >
class HSetWithCones : public HSetT, public QFormT {
public:
  typedef HSetT HSetType;
  typedef QFormT QFormType;  ///< type of quadratic form
  HSetWithCones(const HSetType & hset, const QFormType & Q)
  : HSetType(hset), QFormType(Q){
  }
  HSetWithCones(const HSetType & hset) : HSetType(hset), QFormType(){
  }
  HSetWithCones operator= (const HSetType & hset){
    HSetType::operator=(hset);
    return *this;
  }
  HSetWithCones operator= (const QFormType & Q){
    QFormType::operator=(Q);
    return *this;
  }
  const QFormType & Q() const {
    return (*this);
  }
  QFormType & Q() {
    return (*this);
  }
  virtual std::string show(void) const{
    return std::string("HSet2DWithCones : ") + this->showInfo() + QFormType::show() + "\n";
  }
};
/// @}
}} // end of namespace capd::covrel

#endif /* _CAPD_COVREL_HSETWITHCONES_H_ */
