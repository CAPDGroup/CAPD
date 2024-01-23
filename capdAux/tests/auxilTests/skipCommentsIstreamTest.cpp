/////////////////////////////////////////////////////////////////////////////
/// @file skipCommentsIstream.cpp
///
/// @author Mateusz Juda <mateusz.juda@gmail.com>
///
/// @date 2013-12-01
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD::RedHom Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.edu.pl/ for details.

#include "capd/auxil/skipCommentsIstream.h"

#include <sstream>
#include <iostream>
#include <assert.h>
#include <cstdlib>

#   define ASSERT(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
	      std::exit(1);						\
        } \
    } while (false)

int main()
{
  std::istringstream input("a\n#c1\nb\n%c2\nc\n//\n");
  SkipCommentsIstream stream(input);

  char c;

  stream >> c;
  ASSERT(c == 'a', "c is " << c);
  stream >> c;
  ASSERT(c == 'b', "c is " << c);
  stream >> c;
  ASSERT(c == 'c', "c is " << c);

  return 0;
}
