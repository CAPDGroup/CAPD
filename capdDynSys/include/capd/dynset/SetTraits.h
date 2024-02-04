/////////////////////////////////////////////////////////////////////////////
/// @file SetTraits.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_SETTRAITS_H_
#define _CAPD_DYNSET_SETTRAITS_H_


namespace capd{
namespace dynset{
/// @addtogroup dynset
/// @{

/**
 * This class provides a trait of being set of a given type, i.e. C0Set, C1Set, C2Set and CnSet
 * Used to avoid late binding of virtual function move
 *
*/
template<class SetT>
struct SetTraits{
	const static bool isC0Set=false;
	const static bool isC1Set=false;
	const static bool isC2Set=false;
	const static bool isCnSet=false;
};

/// @}
}} // namespace capd::dynset

#endif

