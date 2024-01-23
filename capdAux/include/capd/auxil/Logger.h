/////////////////////////////////////////////////////////////////////////////
/// @file Logger.h
///
/// @author Mateusz Juda <mateusz.juda@gmail.com>
///
/// @date 2014-06-04
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.edu.pl/ for details.

#ifndef CAPD_FILE_LOGGER_H
#define CAPD_FILE_LOGGER_H



#include <sstream>
#include <iostream>
#include <iterator>
#include <vector>
#include <set>
#include <algorithm>
#include <utility>

namespace capd
{

  namespace auxil
  {

    class Logger {

    public:
      Logger(const std::string& fileName, const std::string& buildDir="", bool global=false);

      std::string name() const { return _name; }

      static void* init()
      {
	static bool init = false;
	if (! init) {
	  Logger::doInit();
	}
	init = true;
	return NULL;
      }

      bool isInfoEnabled() const;
      bool isWarnEnabled() const;
      bool isErrorEnabled() const;
      bool isFatalEnabled() const;
      bool isTraceEnabled() const;
      bool isDebugEnabled() const;

      void forcedLogInfo(const std::string& msg, const char * const fileName,
		     const char * const functionName,
		     int lineNumber);
      void forcedLogWarn(const std::string& msg, const char * const fileName,
		     const char * const functionName,
		     int lineNumber);
      void forcedLogError(const std::string& msg, const char * const fileName,
		     const char * const functionName,
		     int lineNumber);
      void forcedLogFatal(const std::string& msg, const char * const fileName,
		     const char * const functionName,
		     int lineNumber);
      void forcedLogTrace(const std::string& msg, const char * const fileName,
		     const char * const functionName,
		     int lineNumber);
      void forcedLogDebug(const std::string& msg, const char * const fileName,
		     const char * const functionName,
		     int lineNumber);

     private:
      static void doInit();

      std::string _name;
      void* _logger;
    };

  }

}

namespace std {

    template<typename T>
    std::ostream& operator<<(std::ostream& out, const std::vector<T>& vec);

    template<typename T>
    std::ostream& operator<<(std::ostream& out, const std::set<T>& vec);

    template<typename F, typename S>
    std::ostream& operator<<(std::ostream& s, const std::pair<F, S>& p);

    template<typename F, typename S>
    std::ostream& operator<<(std::ostream& s, const std::pair<F, S>& p)
    {
      s << '<' << p.first << ',' << p.second << '>';
      return s;
    }

    namespace impl
    {
      template<typename It>
      std::ostream& printRng(std::ostream& out, It begin, It end)
      {
        out << "{";
        for (It it = begin; it != end; ++it) {
          out << *it << ", ";
        }
        out << "}";

        return out;
      }
    }

    template<typename T>
    std::ostream& operator<<(std::ostream& out, const std::vector<T>& vec)
    {
      return impl::printRng(out, vec.begin(), vec.end());
    }

    template<typename T>
    std::ostream& operator<<(std::ostream& out, const std::set<T>& set)
    {
      return impl::printRng(out, set.begin(), set.end());
    }
}
//   }
// }

#ifndef FILE_BUILD_DIR
   #define FILE_BUILD_DIR ""
//   #error Please define FILE_BUILD_DIR to abs_builddir in Makefile.am (or autoreconf && ./configure)
#endif


#define INIT_CAPD_LOGGER __attribute__((unused)) const void* CAPD_LOGGER_INITIALIZED = capd::auxil::Logger::init();

#define CAPD_CLASS_LOGGER \
  struct CAPD_LOGGER {\
      typedef capd::auxil::Logger& return_type; \
      static inline return_type getCAPDLogger(const char* file, const char* buildDir) { static capd::auxil::Logger logger(file, buildDir, false); return logger; } \
  };

// for global context
struct CAPD_LOGGER {
  typedef capd::auxil::Logger return_type;
  static return_type getCAPDLogger(const char* file, const char* buildDir);
};


// Bellow are macros defined as in LOG4CXX library

#if !defined(CAPD_UNLIKELY)
#if __GNUC__ >= 3
#define CAPD_UNLIKELY(expr) __builtin_expect(expr, 0)
#else
#define CAPD_UNLIKELY(expr) expr
#endif
#endif

  #if !defined(CAPD_LOCATION)
#if defined(_MSC_VER)
#if _MSC_VER >= 1300
      #define __CAPD_FUNC__ __FUNCSIG__
#endif
#else
#if defined(__GNUC__)
      #define __CAPD_FUNC__ __PRETTY_FUNCTION__
#endif
#endif
#if !defined(__CAPD_FUNC__)
#define __CAPD_FUNC__ ""
#endif
      #define CAPD_LOCATION __FILE__, __CAPD_FUNC__, __LINE__
  #endif


#define CAPD_DEBUG(message) { \
        typename CAPD_LOGGER::return_type logger = CAPD_LOGGER::getCAPDLogger(__FILE__, FILE_BUILD_DIR); \
        if (CAPD_UNLIKELY(logger.isDebugEnabled())) {\
           std::stringstream oss_; \
           /*using log4cxx::helpers::operator<<;*/ oss_ << message;     \
           logger.forcedLogDebug(oss_.str(), CAPD_LOCATION); }}

#define CAPD_TRACE(message) { \
        typename CAPD_LOGGER::return_type logger = CAPD_LOGGER::getCAPDLogger(__FILE__, FILE_BUILD_DIR); \
        if (CAPD_UNLIKELY(logger.isTraceEnabled())) {\
           std::stringstream oss_; \
           /*using log4cxx::helpers::operator<<;*/ oss_ << message; \
           logger.forcedLogTrace(oss_.str(), CAPD_LOCATION); }}


#define CAPD_INFO(message) { \
        typename CAPD_LOGGER::return_type logger = CAPD_LOGGER::getCAPDLogger(__FILE__, FILE_BUILD_DIR); \
        if (logger.isInfoEnabled()) {\
           std::stringstream oss_; \
           /*using log4cxx::helpers::operator<<;*/ oss_ << message; \
           logger.forcedLogInfo(oss_.str(), CAPD_LOCATION); }}

#define CAPD_WARN(message) { \
        typename CAPD_LOGGER::return_type logger = CAPD_LOGGER::getCAPDLogger(__FILE__, FILE_BUILD_DIR); \
        if (logger.isWarnEnabled()) {\
           std::stringstream oss_; \
           /*using log4cxx::helpers::operator<<;*/ oss_ << message; \
           logger.forcedLogWarn(oss_.str(), CAPD_LOCATION); }}

#define CAPD_ERROR(message) { \
        typename CAPD_LOGGER::return_type logger = CAPD_LOGGER::getCAPDLogger(__FILE__, FILE_BUILD_DIR); \
        if (logger.isErrorEnabled()) {\
           std::stringstream oss_; \
           /*using log4cxx::helpers::operator<<;*/ oss_ << message; \
           logger.forcedLogError(oss_.str(), CAPD_LOCATION); }}

#define CAPD_FATAL(message) { \
        typename CAPD_LOGGER::return_type logger = CAPD_LOGGER::getCAPDLogger(__FILE__, FILE_BUILD_DIR); \
        if (logger.isFatalEnabled()) {\
           std::stringstream oss_; \
           /*using log4cxx::helpers::operator<<;*/ oss_ << message; \
           logger.forcedLogFatal(oss_.str(), CAPD_LOCATION); }}



// Not required anymore !!!
#define CAPD_DEBUG_TMPL(message) \
       CAPD_DEBUG(message)

#define CAPD_TRACE_TMPL(message) \
       CAPD_TRACE(message)

#define CAPD_INFO_TMPL(message) \
      CAPD_INFO(message)
#define CAPD_WARN_TMPL(message) \
      CAPD_WARN(message)

#define CAPD_ERROR_TMPL(message) \
      CAPD_ERROR(message)

#define CAPD_FATAL_TMPL(message) \
   CAPD_FATAL(message)


#endif // CAPD_FILE_LOGGER_H
