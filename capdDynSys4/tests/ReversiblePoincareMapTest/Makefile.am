include ${capdMake}/make/common_programs_makefile.mkf

check_PROGRAMS= ReversiblePoincareMapTest
TESTS=${check_PROGRAMS}


ReversiblePoincareMapTest_SOURCES=ReversiblePoincareMapTest.cpp
ReversiblePoincareMapTest_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS} ${BOOST_CPPFLAGS}
ReversiblePoincareMapTest_libs=${capdDynSys_LIBS} ${BOOST_UNIT_TEST_FRAMEWORK_LIBS}
ReversiblePoincareMapTest_LDFLAGS=${BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS}
ReversiblePoincareMapTest_LDADD=${ReversiblePoincareMapTest_libs}
ReversiblePoincareMapTest_DEPENDENCIES=${filter %.a %.la,${ReversiblePoincareMapTest_libs}}
