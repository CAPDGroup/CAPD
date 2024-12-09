/////////////////////////////////////////////////////////////////////////////
/// @file powerThree.h
///
/// @author Marian Mrozek
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2006 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.wsb-nlu.edu.pl/ for details.

#ifndef CAPD_POWERTHREE_H
#define CAPD_POWERTHREE_H

#include <iostream>

/// @addtogroup basicalg
/// @{
template<int N>
struct powerThree{
  static const int value=3*powerThree<N-1>::value;
};
/// @}

template<>
struct powerThree<0>{
  static const int value=1;
};

#endif // CAPD_POWERTHREE_H
