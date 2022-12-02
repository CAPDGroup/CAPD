/////////////////////////////////////////////////////////////////////////////
/// @file Logger.cpp
///
/// @author Mateusz Juda <mateusz.juda@gmail.com>
///
/// @date 2014-06-04
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD::RedHom Group.
//
// This file constitutes a part of the CAPD::RedHom library,
// distributed under the terms of the GNU General Public License.
// Consult  http://redhom.ii.edu.pl/ for details.

#include <capd/auxil/Logger.h>
#include <capd/config-capdAux.h>

#if defined(HAVE_BOOST) && defined(HAVE_BOOST_FILESYSTEM) && defined(HAVE_BOOST_SYSTEM) && defined(HAVE_BOOST_REGEX)
#define USE_BOOST 1
#endif

#ifdef USE_BOOST
#include <boost/assign.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/regex.hpp>
#endif

#include <unistd.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

using namespace capd::auxil;

namespace
{
#ifdef USE_BOOST
namespace b = boost;
namespace fs = boost::filesystem;
namespace ba = boost::assign;

std::string getHome()
{
   const char* env = getenv("HOME");
   return env ? env : "/";
}

std::string getConfigFile()
{
#ifndef DATADIR
#error Please define DATADIR as compiler argument
#endif
   fs::path configFile = "capd.log4cxx.properties";

   std::vector<fs::path> paths = ba::list_of(fs::current_path())(fs::path(getHome()))(fs::path(DATADIR))();

   std::stringstream msg;

   msg << "Cannot find " << configFile.string() << " in:";

   BOOST_FOREACH (const fs::path path, paths) {
      if (fs::exists(path / configFile)) {
         return (path / configFile).string();
      } else {
         msg << " " << path.string();
      }
   }

   throw std::runtime_error(msg.str());
}

struct CAPDLoggerName {
   std::string operator()(const std::string& filename_, const std::string& buildDir, bool global)
   {
      // const std::string filename = fs::path(filename).make_preferred().string();

      const bool absolute = fs::path(filename_).is_absolute();
      const std::string filename =
         buildDir.empty() || absolute ? filename_ : (fs::path(buildDir) / fs::path(filename_)).string();

      const std::string defaultName = "capd.default";
      std::string base = "capd";

      const std::string package = this->package(filename);
      const std::string module = this->module(filename, package);
      const std::string tests = this->tests(filename, package);
      const std::string executable = this->executable(filename, package);
      const std::string file = this->file(filename, package);
      const std::string globalSuf = (global ? ".global" : "");

      if (package.empty()) {
         return defaultName;
      } else if (!module.empty()) {
         if (!file.empty()) {
            return base + "." + package + ".lib." + module + "." + file + globalSuf;
         } else {
            return base + "." + package + ".lib." + module + globalSuf;
         }
      } else if (!tests.empty()) {
         if (!file.empty()) {
            return base + "." + package + ".tests." + tests + "." + file + globalSuf;
         } else {
            return base + "." + package + ".tests." + tests + globalSuf;
         }
      } else if (!executable.empty()) {
         if (!file.empty()) {
            return base + "." + package + ".executable." + executable + "." + file + globalSuf;
         } else {
            return base + "." + package + ".executable." + executable + globalSuf;
         }
      }

      return defaultName + globalSuf;
   }

   std::string package(const std::string& filename)
   {
      const b::regex packageRe(".*/(capd(Alg|Aux|Ext|DynSys|RedHom).*?)/.*");
      b::cmatch match;

      if (b::regex_match(filename.c_str(), match, packageRe)) {
         return match[1];
      } else {
         return "";
      }
   }

   std::string module(const std::string& filename, const std::string& package)
   {
      const b::regex moduleRe(".*/" + package + "/(src|include)/(mp)?capd/(.*?)/.*");
      b::cmatch match;

      if (b::regex_match(filename.c_str(), match, moduleRe)) {
         return match[3];
      } else {
         return "";
      }
   }

   std::string tests(const std::string& filename, const std::string& package)
   {
      const b::regex re(".*/" + package + "/tests/(.*?)/.*");
      b::cmatch match;

      if (b::regex_match(filename.c_str(), match, re)) {
         return match[1];
      } else {
         return "";
      }
   }

   std::string executable(const std::string& filename, const std::string& package)
   {
      const b::regex re(".*/" + package + "/(examples|programs)/(.*?)/.*");
      b::cmatch match;

      if (b::regex_match(filename.c_str(), match, re)) {
         return match[2];
      } else {
         return "";
      }
   }

   std::string file(const std::string& filename, const std::string& package)
   {
      const b::regex fileRe(".*/" + package + "/.*/(.*?)\\..*");
      b::cmatch match;

      if (b::regex_match(filename.c_str(), match, fileRe)) {
         return match[1];
      } else {
         return "";
      }
   }
};

#else
std::string getConfigFile() { return ""; }
struct CAPDLoggerName {
   std::string operator()(const std::string& /*filename_*/, const std::string& /*buildDir*/, bool global)
   {
      return std::string("capd.default") + (global ? ".global" : "");
   }
};

#endif  // USE_BOOST
}

#ifdef HAVE_LOG4CXX

#include <log4cxx/basicconfigurator.h>
#include <log4cxx/helpers/exception.h>
#include <log4cxx/logger.h>
#include <log4cxx/propertyconfigurator.h>

#include <memory>

namespace capdLoggerHelper{
    template <typename T>
    T getPointer(T ptr ){
        return ptr;
    }

    template <typename T>
T* getPointer(std::shared_ptr<T> ptr ){
  return ptr.get();
}

}

Logger::Logger(const std::string& fileName, const std::string& buildDir, bool global)
{
   _name = CAPDLoggerName()(fileName, buildDir, global);
   _logger = capdLoggerHelper::getPointer(log4cxx::Logger::getLogger(_name.c_str()));
}

void Logger::doInit()
{
   bool initialized = false;

   try {
      log4cxx::LoggerPtr rootLogger = log4cxx::Logger::getRootLogger();
      initialized = rootLogger->getAllAppenders().size() ? true : false;
      if (initialized) {
         LOG4CXX_DEBUG(rootLogger, "log4cxx initialized using Default Initialization Procedure.");
      }
   } catch (log4cxx::helpers::Exception& e) {
      std::cerr << e.what() << std::endl;
   }

   if (initialized) {
      return;
   }

   std::string configFile;
   std::stringstream msg;
   try {
      configFile = getConfigFile();
   } catch (const std::exception& ex) {
      msg << "\n"
          << ex.what() << "\n"
          << "Please provide above file from CAPD sources (capdAux/src/capd/auxil/) or initialize log4cxx library as "
             "described in the library manual (see http://logging.apache.org/log4cxx/usage.html)."
          << std::endl;
   }

   if (configFile.empty()) {
      log4cxx::BasicConfigurator::configure();
      LOG4CXX_DEBUG(log4cxx::Logger::getRootLogger(),
                    "log4cxx initialized using log4cxx::BasicConfigurator::configure()" << msg.str());
   } else {
      log4cxx::PropertyConfigurator::configure(configFile.c_str());
      LOG4CXX_DEBUG(log4cxx::Logger::getRootLogger(), "log4cxx initialized using: " << configFile);
   }
}

bool Logger::isInfoEnabled() const { return static_cast<log4cxx::Logger*>(_logger)->isInfoEnabled(); }
bool Logger::isWarnEnabled() const { return static_cast<log4cxx::Logger*>(_logger)->isWarnEnabled(); }
bool Logger::isErrorEnabled() const { return static_cast<log4cxx::Logger*>(_logger)->isErrorEnabled(); }
bool Logger::isFatalEnabled() const { return static_cast<log4cxx::Logger*>(_logger)->isFatalEnabled(); }
bool Logger::isTraceEnabled() const { return static_cast<log4cxx::Logger*>(_logger)->isTraceEnabled(); }
bool Logger::isDebugEnabled() const { return static_cast<log4cxx::Logger*>(_logger)->isDebugEnabled(); }
void Logger::forcedLogInfo(const std::string& msg, const char* const fileName, const char* const functionName,
                           int lineNumber)
{
   static_cast<log4cxx::Logger*>(_logger)->forcedLog(::log4cxx::Level::getInfo(), msg,
                                                     ::log4cxx::spi::LocationInfo(fileName, functionName, lineNumber));
}

void Logger::forcedLogWarn(const std::string& msg, const char* const fileName, const char* const functionName,
                           int lineNumber)
{
   static_cast<log4cxx::Logger*>(_logger)->forcedLog(::log4cxx::Level::getWarn(), msg,
                                                     ::log4cxx::spi::LocationInfo(fileName, functionName, lineNumber));
}

void Logger::forcedLogError(const std::string& msg, const char* const fileName, const char* const functionName,
                            int lineNumber)
{
   static_cast<log4cxx::Logger*>(_logger)->forcedLog(::log4cxx::Level::getError(), msg,
                                                     ::log4cxx::spi::LocationInfo(fileName, functionName, lineNumber));
}

void Logger::forcedLogFatal(const std::string& msg, const char* const fileName, const char* const functionName,
                            int lineNumber)
{
   static_cast<log4cxx::Logger*>(_logger)->forcedLog(::log4cxx::Level::getFatal(), msg,
                                                     ::log4cxx::spi::LocationInfo(fileName, functionName, lineNumber));
}

void Logger::forcedLogTrace(const std::string& msg, const char* const fileName, const char* const functionName,
                            int lineNumber)
{
   static_cast<log4cxx::Logger*>(_logger)->forcedLog(::log4cxx::Level::getTrace(), msg,
                                                     ::log4cxx::spi::LocationInfo(fileName, functionName, lineNumber));
}

void Logger::forcedLogDebug(const std::string& msg, const char* const fileName, const char* const functionName,
                            int lineNumber)
{
   static_cast<log4cxx::Logger*>(_logger)->forcedLog(::log4cxx::Level::getDebug(), msg,
                                                     ::log4cxx::spi::LocationInfo(fileName, functionName, lineNumber));
}

#else

Logger::Logger(const std::string& fileName, const std::string& buildDir, bool global)
{
   _name = CAPDLoggerName()(fileName, buildDir, global);
   _logger = &std::cout;
}

void Logger::doInit() {}
bool Logger::isInfoEnabled() const { return true; }
bool Logger::isWarnEnabled() const { return true; }
bool Logger::isErrorEnabled() const { return true; }
bool Logger::isFatalEnabled() const { return true; }
bool Logger::isTraceEnabled() const { return false; }
bool Logger::isDebugEnabled() const { return false; }
void Logger::forcedLogInfo(const std::string& msg, const char* const /*fileName*/, const char* const /*functionName*/,
                           int /*lineNumber*/
                           )
{
   *static_cast<std::ostream*>(_logger) << "INFO " << _name << ": " << msg << std::endl;
}

void Logger::forcedLogWarn(const std::string& msg, const char* const /*fileName*/, const char* const /*functionName*/,
                           int /*lineNumber*/
                           )
{
   *static_cast<std::ostream*>(_logger) << "WARN " << _name << ": " << msg << std::endl;
}

void Logger::forcedLogError(const std::string& msg, const char* const /*fileName*/, const char* const /*functionName*/,
                            int /*lineNumber*/
                            )
{
   *static_cast<std::ostream*>(_logger) << "ERROR " << _name << ": " << msg << std::endl;
}

void Logger::forcedLogFatal(const std::string& msg, const char* const /*fileName*/, const char* const /*functionName*/,
                            int /*lineNumber*/
                            )
{
   *static_cast<std::ostream*>(_logger) << "FATAL " << _name << ": " << msg << std::endl;
}

void Logger::forcedLogTrace(const std::string& msg, const char* const /*fileName*/, const char* const /*functionName*/,
                            int /*lineNumber*/
                            )
{
   *static_cast<std::ostream*>(_logger) << "TRACE " << _name << ": " << msg << std::endl;
}

void Logger::forcedLogDebug(const std::string& msg, const char* const /*fileName*/, const char* const /*functionName*/,
                            int /*lineNumber*/
                            )
{
   *static_cast<std::ostream*>(_logger) << "DEBUG " << _name << ": " << msg << std::endl;
}

#endif  // HAVE_LOG4CXX

// Bellow the logger instance cannot be static. We need different loggers for each file and we have only one global
// function.
// It cost much: we parse file path, build packaga/module/file logger name, hash lookup.
// We do not care too much about it because we should use global looger. In the new code everything should be in a
// class.
capd::auxil::Logger CAPD_LOGGER::getCAPDLogger(const char* file, const char* buildDir)
{
   capd::auxil::Logger logger(file, buildDir, true);
   return logger;
}
