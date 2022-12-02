/////////////////////////////////////////////////////////////////////////////
/// @file mpcapd/diffAlgebra/Hessian.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.


#include "capd/multiPrec/mplib.h"
#include "capd/vectalg/mplib.h"
#include "capd/diffAlgebra/Hessian.hpp"

namespace capd{
namespace diffAlgebra{
/// @addtogroup diffAlgebra 
/// @{
template class capd::diffAlgebra::Hessian<capd::MpFloat,0,0>;
template class capd::diffAlgebra::Hessian<capd::MpInterval,0,0>;

template capd::diffAlgebra::Hessian<capd::MpFloat,0,0> operator*<capd::MpFloat,0>(const capd::MpMatrix&, const capd::diffAlgebra::Hessian<capd::MpFloat,0,0>&);
template capd::diffAlgebra::Hessian<capd::MpInterval,0,0> operator*<capd::MpInterval,0>(const capd::MpIMatrix&, const capd::diffAlgebra::Hessian<capd::MpInterval,0,0>&);

template capd::diffAlgebra::Hessian<capd::MpFloat,0,0> operator*<capd::MpFloat,0>(const capd::diffAlgebra::Hessian<capd::MpFloat,0,0>&, const capd::MpMatrix&);
template capd::diffAlgebra::Hessian<capd::MpInterval,0,0> operator*<capd::MpInterval,0>(const capd::diffAlgebra::Hessian<capd::MpInterval,0,0>&, const capd::MpIMatrix&);

template capd::diffAlgebra::Hessian<capd::MpFloat,0,0> operator*<capd::MpFloat,0>(capd::MpFloat, const capd::diffAlgebra::Hessian<capd::MpFloat,0,0>&);
template capd::diffAlgebra::Hessian<capd::MpInterval,0,0> operator*<capd::MpInterval,0>(capd::MpInterval, const capd::diffAlgebra::Hessian<capd::MpInterval,0,0>&);

}}
