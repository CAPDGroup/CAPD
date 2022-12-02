include ${capdMake}/make/common_programs_makefile.mkf

if HAVE_BOOST
   check_PROGRAMS=C1CurveMaskTest C2CurveMaskTest CnCurveMaskTest
   TESTS=C1CurveMaskTest C2CurveMaskTest CnCurveMaskTest
endif

C1CurveMaskTest_SOURCES=C1CurveMaskTest.cpp compare.h
C1CurveMaskTest_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS} ${BOOST_CPPFLAGS}
C1CurveMaskTest_libs=${capdDynSys_LIBS} ${BOOST_UNIT_TEST_FRAMEWORK_LIBS}
C1CurveMaskTest_LDFLAGS=${BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS}
C1CurveMaskTest_LDADD=${C1CurveMaskTest_libs}
C1CurveMaskTest_DEPENDENCIES=${filter %.a %.la,${C1CurveMaskTest_libs}}

C2CurveMaskTest_SOURCES=C2CurveMaskTest.cpp compare.h
C2CurveMaskTest_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS} ${BOOST_CPPFLAGS}
C2CurveMaskTest_libs=${capdDynSys_LIBS} ${BOOST_UNIT_TEST_FRAMEWORK_LIBS}
C2CurveMaskTest_LDFLAGS=${BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS}
C2CurveMaskTest_LDADD=${C2CurveMaskTest_libs}
C2CurveMaskTest_DEPENDENCIES=${filter %.a %.la,${C2CurveMaskTest_libs}}

CnCurveMaskTest_SOURCES=CnCurveMaskTest.cpp compare.h
CnCurveMaskTest_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS} ${BOOST_CPPFLAGS}
CnCurveMaskTest_libs=${capdDynSys_LIBS} ${BOOST_UNIT_TEST_FRAMEWORK_LIBS}
CnCurveMaskTest_LDFLAGS=${BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS}
CnCurveMaskTest_LDADD=${CnCurveMaskTest_libs}
CnCurveMaskTest_DEPENDENCIES=${filter %.a %.la,${CnCurveMaskTest_libs}}

C1CurveMaskTest_CXXFLAGS+= ${CAPD_CXXFLAGS_CXX11}
C2CurveMaskTest_CXXFLAGS+= ${CAPD_CXXFLAGS_CXX11}
CnCurveMaskTest_CXXFLAGS+= ${CAPD_CXXFLAGS_CXX11}
