include ${capdMake}/make/common_programs_makefile.mkf

check_PROGRAMS=BasicCnTaylorTest
TESTS=${check_PROGRAMS}

BasicCnTaylorTest_SOURCES=BasicCnTaylorTest.cpp
BasicCnTaylorTest_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS} 
BasicCnTaylorTest_libs=${capdDynSys_LIBS} 
#${capdExt_LIBS}
BasicCnTaylorTest_LDADD=${BasicCnTaylorTest_libs}
BasicCnTaylorTest_DEPENDENCIES=${filter %.a %.la,${BasicCnTaylorTest_libs}}

