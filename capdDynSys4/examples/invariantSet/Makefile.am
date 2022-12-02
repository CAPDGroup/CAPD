include ${capdMake}/make/common_programs_makefile.mkf

noinst_PROGRAMS=invSet

invSet_SOURCES=invSet.cpp
invSet_CXXFLAGS= ${AM_CXXFLAGS} ${capdkrak_CFLAGS} ${capdDynSys_CFLAGS} ${capdAux_CFLAGS} -std=c++11 -fopenmp
invSet_libs= ${capdkrak_LIBS} ${capdDynSys_LIBS} ${capdAux_LIBS}  -lpthread
invSet_LDADD=${invSet_libs}
invSet_LDFLAGS=-fopenmp
invSet_DEPENDENCIES=${filter %.a %.la,${invSet_libs}}
