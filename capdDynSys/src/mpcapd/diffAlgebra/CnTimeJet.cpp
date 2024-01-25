/////////////////////////////////////////////////////////////////////////////
/// @file mpcapd/diffAlgebra/CnTimeJet.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.wsb-nlu.edu.pl/ for details.

#include "capd/vectalg/mplib.h"
#include "capd/diffAlgebra/CnTimeJet.h"
#include "capd/diffAlgebra/Jet.hpp"
#include "capd/vectalg/Container.hpp"
#include "capd/vectalg/Matrix.hpp"

template class capd::diffAlgebra::CnContainer<capd::MpFloat,0,0,0>;
template class capd::diffAlgebra::CnContainer<capd::MpInterval,0,0,0>;

template class capd::diffAlgebra::CnContainer<capd::MpFloat,3,3,3>;
template class capd::diffAlgebra::CnContainer<capd::MpInterval,3,3,3>;

template class capd::diffAlgebra::CnTimeJet< capd::MpMatrix, 0 >;
template class capd::diffAlgebra::CnTimeJet< capd::MpIMatrix, 0 >;

template class capd::diffAlgebra::Jet< capd::MpMatrix, 0 >;
template class capd::diffAlgebra::Jet< capd::MpIMatrix, 0 >;

template capd::diffAlgebra::Jet< capd::MpMatrix, 0 > operator+<capd::MpMatrix,0>(const capd::diffAlgebra::Jet< capd::MpMatrix, 0 >&, const capd::diffAlgebra::Jet< capd::MpMatrix, 0 >&);
template capd::diffAlgebra::Jet< capd::MpIMatrix, 0 > operator+<capd::MpIMatrix,0>(const capd::diffAlgebra::Jet< capd::MpIMatrix, 0 >&, const capd::diffAlgebra::Jet< capd::MpIMatrix, 0 >&);

template capd::diffAlgebra::Jet< capd::MpMatrix, 0 > operator-<capd::MpMatrix,0>(const capd::diffAlgebra::Jet< capd::MpMatrix, 0 >&, const capd::diffAlgebra::Jet< capd::MpMatrix, 0 >&);
template capd::diffAlgebra::Jet< capd::MpIMatrix, 0 > operator-<capd::MpIMatrix,0>(const capd::diffAlgebra::Jet< capd::MpIMatrix, 0 >&, const capd::diffAlgebra::Jet< capd::MpIMatrix, 0 >&);

template capd::diffAlgebra::Jet< capd::MpMatrix, 0 > operator*<capd::MpMatrix,0>(const capd::MpMatrix&, const capd::diffAlgebra::Jet< capd::MpMatrix, 0 >&);
template capd::diffAlgebra::Jet< capd::MpIMatrix, 0 > operator*<capd::MpIMatrix,0>(const capd::MpIMatrix&, const capd::diffAlgebra::Jet< capd::MpIMatrix, 0 >&);

template void substitutionPowerSeries<capd::MpMatrix,0>(
      const capd::diffAlgebra::Jet<capd::MpMatrix,0>&,
      const capd::diffAlgebra::Jet<capd::MpMatrix,0>&,
      capd::diffAlgebra::Jet<capd::MpMatrix,0>& result,
      bool
   );
template void substitutionPowerSeries<capd::MpIMatrix,0>(
      const capd::diffAlgebra::Jet<capd::MpIMatrix,0>&,
      const capd::diffAlgebra::Jet<capd::MpIMatrix,0>&,
      capd::diffAlgebra::Jet<capd::MpIMatrix,0>& result,
      bool
   );

template capd::diffAlgebra::Jet<capd::MpMatrix,0> inverseSeriesCloseToIdentity<capd::diffAlgebra::Jet<capd::MpMatrix,0> >(const capd::diffAlgebra::Jet<capd::MpMatrix,0>&);
template capd::diffAlgebra::Jet<capd::MpIMatrix,0> inverseSeriesCloseToIdentity<capd::diffAlgebra::Jet<capd::MpIMatrix,0> >(const capd::diffAlgebra::Jet<capd::MpIMatrix,0>&);

template capd::diffAlgebra::Jet<capd::MpMatrix,0> inversePowerSeries<capd::diffAlgebra::Jet<capd::MpMatrix,0> >(const capd::diffAlgebra::Jet<capd::MpMatrix,0>&, const capd::diffAlgebra::Jet<capd::MpMatrix,0>::MatrixType&);
template capd::diffAlgebra::Jet<capd::MpIMatrix,0> inversePowerSeries<capd::diffAlgebra::Jet<capd::MpIMatrix,0> >(const capd::diffAlgebra::Jet<capd::MpIMatrix,0>&, const capd::diffAlgebra::Jet<capd::MpIMatrix,0>::MatrixType&);
