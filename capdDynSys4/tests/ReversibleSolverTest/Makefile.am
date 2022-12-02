include ${capdMake}/make/common_programs_makefile.mkf

check_PROGRAMS= ReversibleODETest
TESTS=${check_PROGRAMS}


ReversibleODETest_SOURCES = ReversibleODETest.cpp ReversibleSolverTest.cpp FadBadReversibleSolverTest.cpp doTest.h
ReversibleODETest_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS} ${BOOST_CPPFLAGS}
ReversibleODETest_libs=${capdDynSys_LIBS} ${BOOST_UNIT_TEST_FRAMEWORK_LIBS}
ReversibleODETest_LDFLAGS=${BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS}
ReversibleODETest_LDADD=${ReversibleODETest_libs}
ReversibleODETest_DEPENDENCIES=${filter %.a %.la,${ReversibleODETest_libs}}
