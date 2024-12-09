/////////////////////////////////////////////////////////////////////////////
/// @file Vector.h
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.wsb-nlu.edu.pl/ for details.
#ifndef CAPD_VECTALG_EMPTYINTERSECTIONEXCEPTION_H
#define CAPD_VECTALG_EMPTYINTERSECTIONEXCEPTION_H

#include <stdexcept>

/// @addtogroup vectalg
/// @{

class EmptyIntersectionException : public std::runtime_error{
public:
  EmptyIntersectionException(const char * msg):std::runtime_error(msg){
  }
};

/// @}

#endif // CAPD_VECTALG_EMPTYINTERSECTIONEXCEPTION_H
