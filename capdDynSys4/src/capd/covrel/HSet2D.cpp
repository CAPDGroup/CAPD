/// @addtogroup covrel
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file HSet2D.cpp
///
/// @author Daniel Wilczak, Tomasz Kapela
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2016 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include "capd/vectalg/lib.h"
#include "capd/covrel/HSet2D.h"
#include "capd/covrel/GridSet.h"

template class capd::covrel::HSet2D<capd::DMatrix,capd::IMatrix>;

template capd::covrel::GridSet<capd::IMatrix> & capd::covrel::HSet2D<capd::DMatrix,capd::IMatrix>::gridLeftEdge<capd::IMatrix>(
    capd::covrel::GridSet<capd::IMatrix>& G, size_type gridSize, size_type totalDimension = 2, size_type d1=0, size_type d2=1) const;

template capd::covrel::GridSet<capd::IMatrix> & capd::covrel::HSet2D<capd::DMatrix,capd::IMatrix>::gridRightEdge<capd::IMatrix>(
    capd::covrel::GridSet<capd::IMatrix>& G, size_type gridSize, size_type totalDimension = 2, size_type d1=0, size_type d2=1) const;

template capd::covrel::GridSet<capd::IMatrix> & capd::covrel::HSet2D<capd::DMatrix,capd::IMatrix>::gridBottomEdge<capd::IMatrix>(
    capd::covrel::GridSet<capd::IMatrix>& G, size_type gridSize, size_type totalDimension = 2, size_type d1=0, size_type d2=1) const;

template capd::covrel::GridSet<capd::IMatrix> & capd::covrel::HSet2D<capd::DMatrix,capd::IMatrix>::gridTopEdge<capd::IMatrix>(
    capd::covrel::GridSet<capd::IMatrix>& G, size_type gridSize, size_type totalDimension = 2, size_type d1=0, size_type d2=1) const;

template capd::covrel::GridSet<capd::IMatrix> & capd::covrel::HSet2D<capd::DMatrix,capd::IMatrix>::gridSet<capd::IMatrix>(
    capd::covrel::GridSet<capd::IMatrix>& G, size_type grid1, size_type grid2, size_type totalDimension=2, size_type d1=0, size_type d2=1) const;

/// @}
