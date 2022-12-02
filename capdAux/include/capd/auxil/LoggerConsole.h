/////////////////////////////////////////////////////////////////////////////
/// @file LoggerConsole.h
///
/// @author Mateusz Juda <mateusz.juda@gmail.com>
///
/// @date 2014-06-10
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.edu.pl/ for details.

#ifndef CAPD_FILE_LOGGERCONSOLE_H
#define CAPD_FILE_LOGGERCONSOLE_H

#include <capd/config-capdAux.h>

#ifdef HAVE_LOG4CXX
 #error Please use Logger.h
#endif

#include <iostream>

namespace capd
{

  namespace auxil
  {

    class Logger {

    public:
      Logger(const std::string& fileName, const std::string& buildDir="", bool global=false) {}

      std::ostream& operator()()
      {
	return std::cout;
      }
    };

  }

}


#define INIT_CAPD_CXX_LOGGER 0
#define CAPD_CLASS_LOGGER static inline capd::auxil::Logger& getCAPDLogger(const char* file, const char* buildDir) { static capd::auxil::Logger logger(file, buildDir, false); return logger; }

// for global context
capd::auxil::Logger getCAPDLogger(const char* file, const char* buildDir);


#define CAPD_TRACE(message)
#define CAPD_DEBUG(message)
#define CAPD_INFO(message)   ( getCAPDLogger(__FILE__, FILE_BUILD_DIR)() << "INFO: " << message << std::endl)
#define CAPD_WARN(message)   ( getCAPDLogger(__FILE__, FILE_BUILD_DIR)() << "WARN: " <<  message << std::endl)
#define CAPD_ERROR(message) ( getCAPDLogger(__FILE__, FILE_BUILD_DIR)() << "ERROR: " <<  message << std::endl)
#define CAPD_FATAL(message)  ( getCAPDLogger(__FILE__, FILE_BUILD_DIR)() << "FATAL: " << message << std::endl)


#endif // CAPD_FILE_LOGGERCONSOLE_H
