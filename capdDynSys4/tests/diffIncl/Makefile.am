include ${capdMake}/make/common_programs_makefile.mkf

check_PROGRAMS=	diffInclTest
TESTS=${check_PROGRAMS}

diffInclTest_SOURCES=diffInclTest.cpp
diffInclTest_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS} 
diffInclTest_libs=${capdDynSys_LIBS} 
#${capdExt_LIBS}
diffInclTest_LDADD=${diffInclTest_libs}
diffInclTest_DEPENDENCIES=${filter %.a %.la,${diffInclTest_libs}}

