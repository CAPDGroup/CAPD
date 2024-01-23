/////////////////////////////////////////////////////////////////////////////
/// @file Hessian.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.


#include "capd/vectalg/lib.h"
#include "capd/diffAlgebra/Hessian.hpp"

namespace capd{
namespace diffAlgebra{
/// @addtogroup diffAlgebra 
/// @{
template class capd::diffAlgebra::Hessian<double,0,0>;
template class capd::diffAlgebra::Hessian<long double,0,0>;
template class capd::diffAlgebra::Hessian<capd::DInterval,0,0>;

template capd::diffAlgebra::Hessian<double,0,0> operator*<double,0>(const capd::DMatrix&, const capd::diffAlgebra::Hessian<double,0,0>&);
template capd::diffAlgebra::Hessian<long double,0,0> operator*<long double,0>(const capd::LDMatrix&, const capd::diffAlgebra::Hessian<long double,0,0>&);
template capd::diffAlgebra::Hessian<capd::DInterval,0,0> operator*<capd::DInterval,0>(const capd::IMatrix&, const capd::diffAlgebra::Hessian<capd::DInterval,0,0>&);

template capd::diffAlgebra::Hessian<double,0,0> operator*<double,0>(const capd::diffAlgebra::Hessian<double,0,0>&, const capd::DMatrix&);
template capd::diffAlgebra::Hessian<long double,0,0> operator*<long double,0>(const capd::diffAlgebra::Hessian<long double,0,0>&, const capd::LDMatrix&);
template capd::diffAlgebra::Hessian<capd::DInterval,0,0> operator*<capd::DInterval,0>(const capd::diffAlgebra::Hessian<capd::DInterval,0,0>&, const capd::IMatrix&);

template capd::diffAlgebra::Hessian<double,0,0> operator*<double,0>(double, const capd::diffAlgebra::Hessian<double,0,0>&);
template capd::diffAlgebra::Hessian<long double,0,0> operator*<long double,0>(long double, const capd::diffAlgebra::Hessian<long double,0,0>&);
template capd::diffAlgebra::Hessian<capd::DInterval,0,0> operator*<capd::DInterval,0>(capd::DInterval, const capd::diffAlgebra::Hessian<capd::DInterval,0,0>&);

}}
