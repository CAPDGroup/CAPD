/////////////////////////////////////////////////////////////////////////////
/// @file BuildInfoTest
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


#include <capd/auxil/BuildInfo.h>
#include <capd/auxil/Logger.h>

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/typeof/typeof.hpp>

#include <map>
#include <vector>
#include <algorithm>
#include <utility>
#include <sstream>

namespace bl = boost::lambda;
namespace b = boost;
namespace ba = boost::assign;

using namespace capd;
using namespace capd::auxil;

BOOST_AUTO_TEST_SUITE(BuildInfoTestSuite)

BOOST_AUTO_TEST_CASE(output)
{
  std::stringstream out;
  BuildInfo buildInfo;

  out << buildInfo;

  BOOST_CHECK(out.str().size() > 0);
}


BOOST_AUTO_TEST_SUITE_END()
