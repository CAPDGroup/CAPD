/////////////////////////////////////////////////////////////////////////////
/// @file BuildInfo.cpp
///
/// @author Mateusz Juda <mateusz.juda@gmail.com>
///
/// @date 2014-09-23
/////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library (capdAux),
// distributed under the terms of the GNU General Public License.
// Consult http://capd.ii.uj.edu.pl and  http://redhom.ii.edu.pl/ for details.
/////////////////////////////////////////////////////////////////////////////

#include <capd/auxil/BuildInfo.h>
#include <capd/config-capdAux.h>
#include <iostream>

namespace capd
{
  namespace auxil
  {

    const std::string BuildInfo::SVN_REVISION_1 = DEF_SVN_REVISION_1;
    const std::string BuildInfo::SVN_REVISION_2 = DEF_SVN_REVISION_2;
    const std::string BuildInfo::SVN_REVISION_3 = DEF_SVN_REVISION_3;
    const std::string BuildInfo::SVN_REVISION_4 = DEF_SVN_REVISION_4;
    const std::string BuildInfo::SVN_URL_1 = DEF_SVN_URL_1;
    const std::string BuildInfo::SVN_URL_2 = DEF_SVN_URL_2;
    const std::string BuildInfo::SVN_URL_3 = DEF_SVN_URL_3;
    const std::string BuildInfo::SVN_URL_4 = DEF_SVN_URL_4;
    const std::string BuildInfo::BUILD_TAG = DEF_BUILD_TAG;
    const std::string BuildInfo::JOB_NAME = DEF_JOB_NAME;
    const std::string BuildInfo::BUILD_DISPLAY_NAME = DEF_BUILD_DISPLAY_NAME;
    const std::string BuildInfo::BUILD_ID = DEF_BUILD_ID;
    const std::string BuildInfo::NODE_NAME = DEF_NODE_NAME;
    const std::string BuildInfo::BUILD_NUMBER = DEF_BUILD_NUMBER;
    const std::string BuildInfo::VERSION = PACKAGE_VERSION;
    const std::string BuildInfo::BUILD_DATE = DEF_BUILD_DATE;

    #ifdef DEF_LIBRARY_VERSION
    const std::string BuildInfo::LIBRARY_VERSION = DEF_LIBRARY_VERSION;
    #else
    const std::string BuildInfo::LIBRARY_VERSION = "";
    #endif

    std::ostream& operator<<(std::ostream& out, const BuildInfo& buildInfo)
    {
      out << "\n"
        << "SVN_REVISION_1: " << buildInfo.SVN_REVISION_1 << "\n"
        << "SVN_REVISION_2: " << buildInfo.SVN_REVISION_2 << "\n"
        << "SVN_REVISION_3: " << buildInfo.SVN_REVISION_3 << "\n"
        << "SVN_REVISION_4: " << buildInfo.SVN_REVISION_4 << "\n"
        << "SVN_URL_1: " << buildInfo.SVN_URL_1 << "\n"
        << "SVN_URL_2: " << buildInfo.SVN_URL_2 << "\n"
        << "SVN_URL_3: " << buildInfo.SVN_URL_3 << "\n"
        << "SVN_URL_4: " << buildInfo.SVN_URL_4 << "\n"
        << "JOB_NAME: " << buildInfo.JOB_NAME << "\n"
        << "NODE_NAME: " << buildInfo.NODE_NAME << "\n"
        << "BUILD_TAG: " << buildInfo.BUILD_TAG << "\n"
        << "BUILD_DISPLAY_NAME: " << buildInfo.BUILD_DISPLAY_NAME << "\n"
        << "BUILD_ID: " << buildInfo.BUILD_ID << "\n"
        << "BUILD_NUMBER: " << buildInfo.BUILD_NUMBER << "\n"
        << "BUILD_DATE: " << buildInfo.BUILD_DATE << "\n"
        << "VERSION: " << buildInfo.VERSION << "\n"
        << "LIBRARY_VERSION: " << buildInfo.LIBRARY_VERSION << "\n"
        << std::endl;

      return out;
    }

  }
}
