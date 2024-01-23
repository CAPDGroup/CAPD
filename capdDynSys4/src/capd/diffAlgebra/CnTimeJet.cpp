/////////////////////////////////////////////////////////////////////////////
/// @file CnTimeJet.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include "capd/vectalg/lib.h"
#include "capd/diffAlgebra/CnTimeJet.h"
#include "capd/diffAlgebra/Jet.hpp"
#include "capd/vectalg/Container.hpp"
#include "capd/vectalg/ColumnVector.hpp"
#include "capd/vectalg/Matrix.hpp"
#include "capd/vectalg/Vector.hpp"

template class capd::diffAlgebra::CnContainer<double,0,0,0>;
template class capd::diffAlgebra::CnContainer<long double,0,0,0>;
template class capd::diffAlgebra::CnContainer<capd::DInterval,0,0,0>;

template class capd::diffAlgebra::CnContainer<double,3,3,3>;
template class capd::diffAlgebra::CnContainer<long double,3,3,3>;
template class capd::diffAlgebra::CnContainer<capd::DInterval,3,3,3>;

template class capd::diffAlgebra::CnTimeJet< capd::DMatrix, 0 >;
template class capd::diffAlgebra::CnTimeJet< capd::LDMatrix, 0 >;
template class capd::diffAlgebra::CnTimeJet< capd::IMatrix, 0 >;

template class capd::diffAlgebra::Jet< capd::DMatrix, 0 >;
template class capd::diffAlgebra::Jet< capd::LDMatrix, 0 >;
template class capd::diffAlgebra::Jet< capd::IMatrix, 0 >;

template capd::diffAlgebra::Jet< capd::DMatrix, 0 > operator+<capd::DMatrix,0>(const capd::diffAlgebra::Jet< capd::DMatrix, 0 >&, const capd::diffAlgebra::Jet< capd::DMatrix, 0 >&);
template capd::diffAlgebra::Jet< capd::LDMatrix, 0 > operator+<capd::LDMatrix,0>(const capd::diffAlgebra::Jet< capd::LDMatrix, 0 >&, const capd::diffAlgebra::Jet< capd::LDMatrix, 0 >&);
template capd::diffAlgebra::Jet< capd::IMatrix, 0 > operator+<capd::IMatrix,0>(const capd::diffAlgebra::Jet< capd::IMatrix, 0 >&, const capd::diffAlgebra::Jet< capd::IMatrix, 0 >&);

template capd::diffAlgebra::Jet< capd::DMatrix, 0 > operator-<capd::DMatrix,0>(const capd::diffAlgebra::Jet< capd::DMatrix, 0 >&, const capd::diffAlgebra::Jet< capd::DMatrix, 0 >&);
template capd::diffAlgebra::Jet< capd::LDMatrix, 0 > operator-<capd::LDMatrix,0>(const capd::diffAlgebra::Jet< capd::LDMatrix, 0 >&, const capd::diffAlgebra::Jet< capd::LDMatrix, 0 >&);
template capd::diffAlgebra::Jet< capd::IMatrix, 0 > operator-<capd::IMatrix,0>(const capd::diffAlgebra::Jet< capd::IMatrix, 0 >&, const capd::diffAlgebra::Jet< capd::IMatrix, 0 >&);

template capd::diffAlgebra::Jet< capd::DMatrix, 0 > operator*<capd::DMatrix,0>(const capd::DMatrix&, const capd::diffAlgebra::Jet< capd::DMatrix, 0 >&);
template capd::diffAlgebra::Jet< capd::LDMatrix, 0 > operator*<capd::LDMatrix,0>(const capd::LDMatrix&, const capd::diffAlgebra::Jet< capd::LDMatrix, 0 >&);
template capd::diffAlgebra::Jet< capd::IMatrix, 0 > operator*<capd::IMatrix,0>(const capd::IMatrix&, const capd::diffAlgebra::Jet< capd::IMatrix, 0 >&);

template void substitutionPowerSeries<capd::DMatrix,0>(
      const capd::diffAlgebra::Jet<capd::DMatrix,0>&,
      const capd::diffAlgebra::Jet<capd::DMatrix,0>&, 
      capd::diffAlgebra::Jet<capd::DMatrix,0>& result,
      bool
   );
template void substitutionPowerSeries<capd::LDMatrix,0>(
      const capd::diffAlgebra::Jet<capd::LDMatrix,0>&,
      const capd::diffAlgebra::Jet<capd::LDMatrix,0>&,
      capd::diffAlgebra::Jet<capd::LDMatrix,0>& result,
      bool
   );
template void substitutionPowerSeries<capd::IMatrix,0>(
      const capd::diffAlgebra::Jet<capd::IMatrix,0>&,
      const capd::diffAlgebra::Jet<capd::IMatrix,0>&,
      capd::diffAlgebra::Jet<capd::IMatrix,0>& result,
      bool
   );

template capd::diffAlgebra::Jet<capd::DMatrix,0> inverseSeriesCloseToIdentity<capd::diffAlgebra::Jet<capd::DMatrix,0> >(const capd::diffAlgebra::Jet<capd::DMatrix,0>&);
template capd::diffAlgebra::Jet<capd::LDMatrix,0> inverseSeriesCloseToIdentity<capd::diffAlgebra::Jet<capd::LDMatrix,0> >(const capd::diffAlgebra::Jet<capd::LDMatrix,0>&);
template capd::diffAlgebra::Jet<capd::IMatrix,0> inverseSeriesCloseToIdentity<capd::diffAlgebra::Jet<capd::IMatrix,0> >(const capd::diffAlgebra::Jet<capd::IMatrix,0>&);

template capd::diffAlgebra::Jet<capd::DMatrix,0> inversePowerSeries<capd::diffAlgebra::Jet<capd::DMatrix,0> >(const capd::diffAlgebra::Jet<capd::DMatrix,0>&, const capd::diffAlgebra::Jet<capd::DMatrix,0>::MatrixType&);
template capd::diffAlgebra::Jet<capd::LDMatrix,0> inversePowerSeries<capd::diffAlgebra::Jet<capd::LDMatrix,0> >(const capd::diffAlgebra::Jet<capd::LDMatrix,0>&, const capd::diffAlgebra::Jet<capd::LDMatrix,0>::MatrixType&);
template capd::diffAlgebra::Jet<capd::IMatrix,0> inversePowerSeries<capd::diffAlgebra::Jet<capd::IMatrix,0> >(const capd::diffAlgebra::Jet<capd::IMatrix,0>&, const capd::diffAlgebra::Jet<capd::IMatrix,0>::MatrixType&);
