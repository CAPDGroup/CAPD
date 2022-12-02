/////////////////////////////////////////////////////////////////////////////
/// @file logger.cpp
///
/// @author Mateusz Juda <mateusz.juda@gmail.com>
///
/// @date 2014-06-10
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD::RedHom Group.
//
// This file constitutes a part of the CAPD::RedHom library,
// distributed under the terms of the GNU General Public License.
// Consult  http://redhom.ii.edu.pl/ for details.


#include "Example.h"
#include "capd/auxil/Logger.h"


void withoutLog4CXXExample();

void defaultLoggerExample()
{
  Example example;
  Example::logStatic();
  example.log();

  globalFunction();
}


int main(int /*argc*/, char */*argv*/[])
{

  INIT_CAPD_LOGGER;

  defaultLoggerExample();
  return 0;
}
