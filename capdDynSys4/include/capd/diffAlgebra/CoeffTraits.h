/////////////////////////////////////////////////////////////////////////////
/// @file CoeffTraits.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DIFFALGEBRA_COEFFTRAITS_H_
#define _CAPD_DIFFALGEBRA_COEFFTRAITS_H_


namespace capd{
namespace diffAlgebra{
/// @addtogroup diffAlgebra
/// @{

/**
 * This class provides a trait of being set of a given type, i.e. C0Set, C1Set, C2Set and CnSet
 * Used to avoid late binding of virtual function move
 *
*/
template<class CoeffT>
struct CoeffTraits{
	const static bool isC0Jet=false;
	const static bool isC1Jet=false;
	const static bool isC2Jet=false;
	const static bool isCnJet=false;
};
/// @}
}}

#endif



