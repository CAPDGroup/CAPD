include ${capdMake}/make/common_programs_makefile.mkf

if HAVE_BOOST
   check_PROGRAMS=C1MaskTest C2MaskTest CnMaskTest
   TESTS=C1MaskTest C2MaskTest CnMaskTest
endif

C1MaskTest_SOURCES=C1MaskTest.cpp compare.h
C1MaskTest_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS} ${BOOST_CPPFLAGS}
C1MaskTest_libs=${capdDynSys_LIBS} ${BOOST_UNIT_TEST_FRAMEWORK_LIBS}
C1MaskTest_LDFLAGS=${BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS}
C1MaskTest_LDADD=${C1MaskTest_libs}
C1MaskTest_DEPENDENCIES=${filter %.a %.la,${C1MaskTest_libs}}

C2MaskTest_SOURCES=C2MaskTest.cpp compare.h
C2MaskTest_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS} ${BOOST_CPPFLAGS}
C2MaskTest_libs=${capdDynSys_LIBS} ${BOOST_UNIT_TEST_FRAMEWORK_LIBS}
C2MaskTest_LDFLAGS=${BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS}
C2MaskTest_LDADD=${C2MaskTest_libs}
C2MaskTest_DEPENDENCIES=${filter %.a %.la,${C2MaskTest_libs}}

CnMaskTest_SOURCES=CnMaskTest.cpp compare.h
CnMaskTest_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS} ${BOOST_CPPFLAGS}
CnMaskTest_libs=${capdDynSys_LIBS} ${BOOST_UNIT_TEST_FRAMEWORK_LIBS}
CnMaskTest_LDFLAGS=${BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS}
CnMaskTest_LDADD=${CnMaskTest_libs}
CnMaskTest_DEPENDENCIES=${filter %.a %.la,${CnMaskTest_libs}}


    C1MaskTest_CXXFLAGS+= ${CAPD_CXXFLAGS_CXX11}
    C2MaskTest_CXXFLAGS+= ${CAPD_CXXFLAGS_CXX11}
    CnMaskTest_CXXFLAGS+= ${CAPD_CXXFLAGS_CXX11}

