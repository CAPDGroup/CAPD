include ${capdMake}/make/common_programs_makefile.mkf

if HAVE_BOOST
   check_PROGRAMS=DPoincareMapMaskTest
   TESTS=DPoincareMapMaskTest
endif

DPoincareMapMaskTest_SOURCES=DPoincareMapMaskTest.cpp 
DPoincareMapMaskTest_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS} ${BOOST_CPPFLAGS}
DPoincareMapMaskTest_libs=${capdDynSys_LIBS} ${BOOST_UNIT_TEST_FRAMEWORK_LIBS}
DPoincareMapMaskTest_LDFLAGS=${BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS}
DPoincareMapMaskTest_LDADD=${DPoincareMapMaskTest_libs}
DPoincareMapMaskTest_DEPENDENCIES=${filter %.a %.la,${DPoincareMapMaskTest_libs}}

    DPoincareMapMaskTest_CXXFLAGS+= ${CAPD_CXXFLAGS_CXX11}
