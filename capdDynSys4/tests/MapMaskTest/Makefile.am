include ${capdMake}/make/common_programs_makefile.mkf

if HAVE_BOOST
   check_PROGRAMS=MapMaskTest
   TESTS=MapMaskTest
endif

MapMaskTest_SOURCES=MapMaskTest.cpp EmptyMaskTest.cpp SetMaskTest.cpp AddMaskTest.cpp
MapMaskTest_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS} ${BOOST_CPPFLAGS}
MapMaskTest_libs=${capdDynSys_LIBS} ${BOOST_UNIT_TEST_FRAMEWORK_LIBS}
MapMaskTest_LDFLAGS=${BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS}
MapMaskTest_LDADD=${MapMaskTest_libs}
MapMaskTest_DEPENDENCIES=${filter %.a %.la,${MapMaskTest_libs}}


    MapMaskTest_CXXFLAGS+= ${CAPD_CXXFLAGS_CXX11}

