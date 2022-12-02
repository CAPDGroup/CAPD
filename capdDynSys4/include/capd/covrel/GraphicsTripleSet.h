/// @addtogroup covrel
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file GraphicsTripleSet.h
///
/// @author Jaroslaw Dulak, Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2016 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_COVREL_GRAPHICSTRIPLESET_H_
#define _CAPD_COVREL_GRAPHICSTRIPLESET_H_

#include "capd/krak/IFrame.h"
#include "capd/covrel/TripleSet.h"


namespace capd{
namespace covrel{

class GraphicsTripleSet : public TripleSet
{
public:
  typedef TripleSet::FloatVector     FloatVector;
  typedef TripleSet::IntervalVector  IntervalVector;

  explicit GraphicsTripleSet(const TripleSet& t);
  GraphicsTripleSet(const FloatVector& a_center, const FloatVector& a_stable, const FloatVector& a_unstable)
    : TripleSet(a_center,a_stable,a_unstable){}
  GraphicsTripleSet(const IntervalVector& a_center, const IntervalVector& a_stable, const IntervalVector& a_unstable)
    : TripleSet(a_center,a_stable,a_unstable){}
  GraphicsTripleSet(){}
  ~GraphicsTripleSet(){}

  friend krak::IFrame& operator << (krak::IFrame&, GraphicsTripleSet&);
  void show(krak::IFrame &f, const FloatVector& point, int color=BLACK) const;
  void show(krak::IFrame &f, const IntervalVector& point, int color=BLACK) const;
  void show(krak::IFrame &f) const;
}; //end of class GraphicsTripleSet

}} //namespace capd::covrel

#endif // _CAPD_COVREL_GRAPHICSTRIPLESET_H_

/// @}
