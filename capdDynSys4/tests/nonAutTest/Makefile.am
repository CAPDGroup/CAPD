include ${capdMake}/make/common_programs_makefile.mkf

check_PROGRAMS=nonAutTest
TESTS=${check_PROGRAMS}

nonAutTest_SOURCES=nonAutTest.cpp
nonAutTest_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS}
nonAutTest_libs=${capdDynSys_LIBS}
nonAutTest_LDADD=${nonAutTest_libs}
nonAutTest_DEPENDENCIES=${filter %.a %.la,${nonAutTest_libs}}

