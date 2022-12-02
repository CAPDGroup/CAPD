include ${capdMake}/make/common_programs_makefile.mkf

check_PROGRAMS= runTests
TESTS=${check_PROGRAMS}


runTests_SOURCES=runTests.cpp HighOrderTest.cpp
runTests_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS} ${BOOST_CPPFLAGS}
runTests_libs=${capdDynSys_LIBS} ${BOOST_UNIT_TEST_FRAMEWORK_LIBS}
runTests_LDFLAGS=${BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS}
runTests_LDADD=${runTests_libs}
runTests_DEPENDENCIES=${filter %.a %.la,${runTests_libs}}
