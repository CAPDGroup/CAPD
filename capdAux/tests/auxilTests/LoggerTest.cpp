/////////////////////////////////////////////////////////////////////////////
/// @file loggerTest.cpp
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

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>

#include <map>
#include <vector>
#include <algorithm>
#include <utility>
#include <utility>

namespace bl = boost::lambda;
namespace ba = boost::assign;
namespace b = boost;

using namespace boost::assign;
using namespace boost::lambda;
using namespace capd::auxil;


BOOST_AUTO_TEST_SUITE(loggerTestSuite)

BOOST_AUTO_TEST_CASE(empty)
{
}

#ifdef HAVE_LOG4CXX

BOOST_AUTO_TEST_CASE( checkLoggerName )
{
  namespace b = boost;
  typedef std::pair<std::string, std::string> FileToLogger;
  std::vector<FileToLogger> fixtures = ba::list_of<FileToLogger>
    ("/home/me/work/capdAux/src/capd/auxil/subdir/logger.cpp", "capd.capdAux.lib.auxil.logger")
    ("/home/me/work/capdAux/src/mpcapd/auxil/subdir/logger.cpp", "capd.capdAux.lib.auxil.logger")
    ("/home/me/work/capdAux/src/capd/auxil/logger.cpp", "capd.capdAux.lib.auxil.logger")
    ("/home/me/work/capdAux/include/capd/auxil/logger.h", "capd.capdAux.lib.auxil.logger")
    ("/home/me/work/capdAux/include/capd/auxil/Logger.hpp", "capd.capdAux.lib.auxil.Logger")
    ("/home/me/work/capdAux/include/capd/auxil/logger.X.Y", "capd.capdAux.lib.auxil.logger")
    ("/home/me/work/capdAux/include/capd/auxil/", "capd.capdAux.lib.auxil")
    (std::string(FILE_BUILD_DIR) + "/" + __FILE__ , "capd.capdAux.tests.auxilTests.LoggerTest")
    ("/home/me/work/capdA/include/capd/auxil/", "capd.default")
    ("/home/me/work/capdAux/include/a/auxil/", "capd.default")
    ("/home/me/work/capdAux/inc/capd/auxil/", "capd.default")
    ("/home/me/work/capdAux/capd/auxil/", "capd.default")

    //    ("c:\\home\\me\\work\\capdAux\\src\\capd\\auxil\\logger.cpp", "capd.capdAux.auxil")
    ;

  BOOST_FOREACH(const FileToLogger& fileToLogger, fixtures) {
    Logger logger(fileToLogger.first);

    //    BOOST_TEST_MESSAGE(fileToLogger.first);
    BOOST_CHECK_EQUAL(fileToLogger.second, logger.name());
  }
}


struct TestStruct
{

  TestStruct()
  {
    CAPD_TRACE("X");
  }

  void f() {}

  CAPD_CLASS_LOGGER;
};

template<typename X>
struct TestStructTemplate
{

  TestStructTemplate()
  {
    CAPD_TRACE_TMPL("X");
  }

  void f() {}

  CAPD_CLASS_LOGGER;
};

BOOST_AUTO_TEST_CASE(LoggerInTemplate) // check if compile
{
  TestStructTemplate<int> x;
  TestStruct y;
  x.f();
  y.f();
}

BOOST_AUTO_TEST_CASE(log)
{
  std::vector<int> vec(10);
  std::pair<int, bool> p;

  CAPD_INFO(vec << " " << vec << p);
}

#endif

BOOST_AUTO_TEST_SUITE_END()
