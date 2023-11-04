/////////////////////////////////////////////////////////////////////////////
/// @file CAPDConfig.h
///
/// @author Mateusz Juda <mateusz.juda@gmail.com>
///
/// @date 2014-06-18
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.edu.pl/ for details.

#ifndef CAPD_FILE_CAPDCONFIG_H
#define CAPD_FILE_CAPDCONFIG_H

#include <memory>

namespace capd
{

  namespace auxil
  {
    class ConfigFileReader;

    class CAPDConfig
    {
      CAPDConfig();

    public:
      static CAPDConfig& getInstance();

      bool usePARI();

    private:
      std::unique_ptr<ConfigFileReader> _configFileReader;
    };
  }

}

#endif // CAPD_FILE_CAPDCONFIG_H
