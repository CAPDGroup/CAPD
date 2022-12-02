include ${capdMake}/make/common_programs_makefile.mkf

check_PROGRAMS=
if HAVE_BOOST
check_PROGRAMS+=IComplexTests DComplexTests
endif

TESTS=${check_PROGRAMS}

IComplexTests_SOURCES=IComplexTests.cpp
IComplexTests_libs=${capdAlg_LIBS} ${capdExt_LIBS}  ${BOOST_UNIT_TEST_FRAMEWORK_LIBS}
IComplexTests_LDADD=${IComplexTests_libs}
IComplexTests_DEPENDENCIES=${filter %.a %.la,${IComplexTests_libs}}
IComplexTests_CXXFLAGS= ${AM_CXXFLAGS} ${capdintervals_CFLAGS} ${capdExt_CFLAGS} ${BOOST_CPPFLAGS}


DComplexTests_SOURCES=DComplexTests.cpp
DComplexTests_libs=${capdAlg_LIBS} ${capdExt_LIBS} ${BOOST_UNIT_TEST_FRAMEWORK_LIBS}
DComplexTests_LDADD=${DComplexTests_libs}
DComplexTests_DEPENDENCIES=${filter %.a %.la,${DComplexTests_libs}}
DComplexTests_CXXFLAGS= ${AM_CXXFLAGS} ${capdintervals_CFLAGS} ${capdExt_CFLAGS} ${BOOST_CPPFLAGS}
