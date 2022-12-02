/// @addtogroup covrel
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file GraphicsHSet3D.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2016 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_COVREL_GRAPHICSHSET3D_H_
#define _CAPD_COVREL_GRAPHICSHSET3D_H_

#include "capd/krak/IFrame.h"
#include "capd/covrel/HSet3D.h"



namespace capd{

  namespace krak {
    class Coord3D;
  }

namespace covrel{

class GraphicsHSet3D : public HSet3D{
public:
  typedef HSet3D::FloatVector FloatVector;
  typedef HSet3D::IntervalVector IntervalVector;
  typedef HSet3D::FloatMatrix FloatMatrix;
  typedef HSet3D::IntervalMatrix IntervalMatrix;

  GraphicsHSet3D(){}
  GraphicsHSet3D(const FloatVector &Center, const FloatVector &U, const FloatVector &S1, const FloatVector &S2)
    : HSet3D(Center,U,S1,S2){}
  GraphicsHSet3D(const IntervalVector &Center, const IntervalVector &U, const IntervalVector &S1, const IntervalVector &S2)
    : HSet3D(Center,U,S1,S2){}
  ~GraphicsHSet3D(){}

  void show(krak::IFrame &f, int dim1, int dim2) const;
  void show(krak::IFrame &f, int dim1, int dim2, const IntervalVector &, int color=BLACK) const;
  void show(krak::Coord3D &f) const;
  void show(krak::Coord3D &f, const IntervalVector &iv, int color=BLACK) const;
  friend krak::Coord3D& operator << (krak::Coord3D &, const GraphicsHSet3D&);
};

}} // namespace capd::covrel

#endif // _CAPD_COVREL_GRAPHICSHSET3D_H_

/// @}
