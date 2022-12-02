include ${capdMake}/make/common_programs_makefile.mkf

check_PROGRAMS= SolverGetStepTest
TESTS=${check_PROGRAMS}


SolverGetStepTest_SOURCES=SolverGetStepTest.cpp DFixStep.cpp DVarStep.cpp IFixStep.cpp IVarStep.cpp checkGetStep.h
SolverGetStepTest_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS} ${BOOST_CPPFLAGS}
SolverGetStepTest_libs=${capdDynSys_LIBS} ${BOOST_UNIT_TEST_FRAMEWORK_LIBS}
SolverGetStepTest_LDFLAGS=${BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS}
SolverGetStepTest_LDADD=${SolverGetStepTest_libs}
SolverGetStepTest_DEPENDENCIES=${filter %.a %.la,${SolverGetStepTest_libs}}
