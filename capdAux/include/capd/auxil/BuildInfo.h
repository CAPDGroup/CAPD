/////////////////////////////////////////////////////////////////////////////
/// @file BuildInfo.h
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

#ifndef CAPD_FILE_CAPDAUX_AUXIL_BUILDINFO_H
#define CAPD_FILE_CAPDAUX_AUXIL_BUILDINFO_H

// #include
#include <string>
#include <iosfwd>

namespace capd
{
  namespace auxil
  {

    class BuildInfo
    {

    public:
      static const std::string SVN_REVISION_1;
      static const std::string SVN_REVISION_2;
      static const std::string SVN_REVISION_3;
      static const std::string SVN_REVISION_4;
      static const std::string SVN_URL_1;
      static const std::string SVN_URL_2;
      static const std::string SVN_URL_3;
      static const std::string SVN_URL_4;
      static const std::string BUILD_TAG;
      static const std::string JOB_NAME;
      static const std::string BUILD_DISPLAY_NAME;
      static const std::string BUILD_ID;
      static const std::string NODE_NAME;
      static const std::string BUILD_NUMBER;
      static const std::string VERSION;
      static const std::string BUILD_DATE;
      static const std::string LIBRARY_VERSION;


    private:
      friend std::ostream& operator<<(std::ostream& out, const BuildInfo& buildInfo);

    };


  }
}

#endif // CAPD_FILE_CAPDAUX_AUXIL_BUILDINFO_H
