/// @addtogroup auxil
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file stringOstream.h
///
/// @author Marian Mrozek
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2006 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.wsb-nlu.edu.pl/ for details. 

#ifndef CAPD_STRINGOSTREAM_H
#define CAPD_STRINGOSTREAM_H

#include <sstream>

template<typename T>
std::string& operator<<(std::string& s,const T& t){
  std::ostringstream is;
  is << t;
  s+=is.str();
  return s;
}

#endif // CAPD_STRINGOSTREAM_H

/// @}
