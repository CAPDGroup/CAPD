/////////////////////////////////////////////////////////////////////////////
/// @file CAPDConfig.cpp
///
/// @author Mateusz Juda <mateusz.juda@gmail.com>
///
/// @date 2014-06-18
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD::RedHom Group.
//
// This file constitutes a part of the CAPD::RedHom library,
// distributed under the terms of the GNU General Public License.
// Consult  http://redhom.ii.edu.pl/ for details.

#include <capd/auxil/CAPDConfig.h>
#include <capd/auxil/ConfigFileReader.h>

#include <fstream>
#include <unistd.h>
#include <string>
#include <pwd.h>

using namespace capd::auxil;

CAPDConfig::CAPDConfig():
  _configFileReader(new ConfigFileReader())
{

  std::string filename = "capd.ini";
  if (access(filename.c_str(), R_OK) == -1) {
    const int myuid = getuid();
    passwd *mypasswd = getpwuid(myuid);
    filename = mypasswd->pw_dir + std::string("/") + filename;
  }

  if (access(filename.c_str(), R_OK) != -1) {
    std::ifstream file(filename.c_str());
    if (file) {
      _configFileReader->parse(file);
    }
  }
}


CAPDConfig& CAPDConfig::getInstance()
{
  static CAPDConfig instance;
  return instance;
}

bool CAPDConfig::usePARI()
{
  return _configFileReader->get<bool>("capdAlg.usePARI", true);
}
