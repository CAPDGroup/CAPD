
/////////////////////////////////////////////////////////////////////////////
/// @file C2DoubletonSet.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include "capd/vectalg/lib.h"
#include "capd/dynset/C2DoubletonSet.hpp"
#include "capd/dynset/reorganization/FactorReorganization.h"
#include "capd/dynset/reorganization/QRReorganization.h"
#include "capd/dynset/QRPolicy.h"

using namespace capd::dynset;

typedef QRReorganization<InverseQRPolicy<> > Pped2Policies;
typedef FactorReorganization<FullQRWithPivoting<> > Rect2Policies;

namespace capd{
namespace dynset{

template class C2DoubletonSet< capd::IMatrix , Pped2Policies >;
template class capd::dynset::C2DoubletonSet< capd::IMatrix , Rect2Policies >;

}}
