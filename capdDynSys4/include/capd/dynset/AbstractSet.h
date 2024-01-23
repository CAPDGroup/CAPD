/////////////////////////////////////////////////////////////////////////////
/// @file AbstractSet.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_ABSTRACTSET_H_
#define _CAPD_DYNSET_ABSTRACTSET_H_

#include <string>

namespace capd{
namespace dynset{
/// @addtogroup dynset 
/// @{
template<typename VectorT>
class AbstractSet
{
public:
  using VectorType = VectorT;
  using ScalarType = typename VectorType::ScalarType;

  /// This method computes value of functor f at interval vector represented by this set.
  template<class Functional>
  ScalarType evalAt(const Functional& f) const {
    return f((VectorType)(*this));
  }

  virtual operator VectorType() const = 0;

  /// destructor
  virtual ~AbstractSet (void) = default;
  /// returns a set detailed information
  virtual std::string show() const = 0;
};

/// @}
}} //namespace capd::dynset

#endif // _CAPD_DYNSET_ABSTRACTSET_H_
