include ${capdMake}/make/common_programs_makefile.mkf

check_PROGRAMS= mapTest
#check_PROGRAMS= acosTest asinTest atanTest cosTest expTest logTest sinTest
TESTS=${check_PROGRAMS}


mapTest_SOURCES=mapTest.cpp acosTest.cpp asinTest.cpp atanTest.cpp cosTest.cpp expTest.cpp logTest.cpp sinTest.cpp sqrTest.cpp cubeTest.cpp quarticTest.cpp natPowTest.cpp compare.h
mapTest_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS} ${BOOST_CPPFLAGS}
mapTest_libs=${capdDynSys_LIBS} ${BOOST_UNIT_TEST_FRAMEWORK_LIBS}
#${capdExt_LIBS}
mapTest_LDFLAGS=${BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS}
mapTest_LDADD=${mapTest_libs}
mapTest_DEPENDENCIES=${filter %.a %.la,${mapTest_libs}}
