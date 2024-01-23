/////////////////////////////////////////////////////////////////////////////
/// @file BuildInfoMain.cpp
///
/// @author Mateusz Juda <mateusz.juda@gmail.com>
///
/// @date 2014-09-25
/////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library (capdAux),
// distributed under the terms of the GNU General Public License.
// Consult http://capd.ii.uj.edu.pl and  http://redhom.ii.edu.pl/ for details.
/////////////////////////////////////////////////////////////////////////////

#include <capd/auxil/ApplicationDesc.h>
#include <capd/auxil/Logger.h>
#include <capd/auxil/BuildInfo.h>

#include <string>
#include <iostream>


using namespace capd;

namespace {
  const capd::auxil::ApplicationDesc applicationDesc("capdAux::example::BuildInfoMain",
                                                     "Mateusz Juda <mateusz.juda@gmail.com>",
                                                     "2014-09-25",
                                                     "");
}

int main(int /*argc*/, char */*argv*/[])
{

  try {
    INIT_CAPD_LOGGER;
    CAPD_INFO("Entering application:\n" << applicationDesc);

    capd::auxil::BuildInfo buildInfo;

    CAPD_INFO(buildInfo);
    return 0;

  } catch (std::exception& ex) {
    std::cerr << "Error: " << ex.what() << std::endl;
    return 1;
  } catch (const char* ex) {
    std::cerr << "Error: " << ex << std::endl;
    return 1;
  } catch (const std::string& ex) {
    std::cerr << "Error: " << ex << std::endl;
    return 1;
  } catch(...) {
    std::cerr << "Unknown exception" << std::endl;
    return 1;
  }

  return 1;
}
