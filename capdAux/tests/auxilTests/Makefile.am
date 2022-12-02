include ${capdMake}/make/common_programs_makefile.mkf

check_PROGRAMS=CRefTest  functorTest skipCommentsIstreamTest

if HAVE_BOOST
  check_PROGRAMS+=auxilTest
endif


TESTS=${check_PROGRAMS}

CRefTest_SOURCES=CRefTest.cpp
CRefTest_CXXFLAGS= ${AM_CXXFLAGS} ${capdauxil_CFLAGS}

functorTest_SOURCES=functorTest.cpp
functorTest_CXXFLAGS= ${AM_CXXFLAGS} ${capdauxil_CFLAGS}

skipCommentsIstreamTest_SOURCES=skipCommentsIstreamTest.cpp
skipCommentsIstreamTest_CXXFLAGS= ${AM_CXXFLAGS} ${capdauxil_CFLAGS}
skipCommentsIstreamTest_LDFLAGS=${capdauxil_LIBS}



auxilTest_SOURCES=auxilTestMain.cpp LoggerTest.cpp BuildInfoTest.cpp
auxilTest_CXXFLAGS= ${AM_CXXFLAGS} ${capdAux_CFLAGS} -I${top_builddir}/../include
auxilTest_libs=${capdAux_LIBS} ${BOOST_UNIT_TEST_FRAMEWORK_LIBS}
auxilTest_LDFLAGS=${BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS}
auxilTest_LDADD=${auxilTest_libs}
auxilTest_DEPENDENCIES=${filter %.a %.la,${auxilTest_libs}}
