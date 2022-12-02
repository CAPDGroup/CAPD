/// @addtogroup covrel
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file GraphicsTripleSet.cpp
///
/// @author Jaroslaw Dulak, Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2016 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include "capd/covrel/GraphicsTripleSet.h"

namespace capd{
namespace covrel{

GraphicsTripleSet::GraphicsTripleSet(const TripleSet& t)
  : TripleSet(t)
{}

void GraphicsTripleSet::show(krak::IFrame &f) const
{
  // plotting in Yellow - the N^- edges - but for better visualition
  FloatVector s2 = 2.*this->stable;
  FloatVector E,S;
  E = S = this->center + this->unstable;
  E += s2;
  S -= s2;
  f.line(S[0],S[1],E[0],E[1],YELLOW);

  E = S = this->center - this->unstable;
  S -= s2;
  E += s2;
  f.line(S[0],S[1],E[0],E[1],YELLOW);

  // plotting the support - parallelogram
  f.jump(A[0],A[1]);
  f.draw(B[0],B[1]);
  f.draw(C[0],C[1]);
  f.draw(D[0],D[1]);
  f.draw(A[0],A[1]);

  // plotting the stable direction in BLUE
  f.jump(this->center[0],this->center[1]);
  f.draw(this->center[0] + this->stable[0],this->center[1] + this->stable[1],BLUE);

  // plotting the unstable direction in RED
  f.jump(this->center[0],this->center[1]);
  f.draw(this->center[0] + this->unstable[0],center[1] + this->unstable[1],RED);
}

// -----------------------------------------------------------------

void GraphicsTripleSet::show(krak::IFrame &f, const FloatVector& point, int color) const
{
  show(f);													//we plot support
  f.Xcrss(point[0],point[1],5,color);  //we plot the point
}

// -----------------------------------------------------------------

void GraphicsTripleSet::show(krak::IFrame &f, const IntervalVector& vect, int color) const
{
  show(f);
  f.boxFill(vect[0].leftBound(),vect[1].leftBound(),vect[0].rightBound(),vect[1].rightBound(),color);
  f.box(vect[0].leftBound(),vect[1].leftBound(),vect[0].rightBound(),vect[1].rightBound(),BLACK);
}

// -----------------------------------------------------------------

krak::IFrame& operator<<(krak::IFrame& f, GraphicsTripleSet& Set)
{
  Set.show(f);
  return f;
}

}} // namespace capd::covrel

/// @}
