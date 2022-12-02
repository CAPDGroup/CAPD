include ${capdMake}/make/common_programs_makefile.mkf

noinst_PROGRAMS=diffInclExample

diffInclExample_SOURCES=diffInclExample.cpp
diffInclExample_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS}
diffInclExample_libs=${capdDynSys_LIBS}
diffInclExample_LDADD=${diffInclExample_libs}
diffInclExample_DEPENDENCIES=${filter %.a %.la,${diffInclExample_libs}}

EXTRA_DIST=diffInclExample.dox
