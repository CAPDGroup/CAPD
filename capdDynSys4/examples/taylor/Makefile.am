include ${capdMake}/make/common_programs_makefile.mkf

noinst_PROGRAMS=tayltst

tayltst_SOURCES=tayltst.cpp
tayltst_CXXFLAGS= ${AM_CXXFLAGS} ${capdkrak_CFLAGS} ${capdDynSys_CFLAGS} ${capdAux_CFLAGS}
tayltst_libs=${capdkrak_LIBS} ${capdDynSys_LIBS} ${capdAux_LIBS}
tayltst_LDADD=${tayltst_libs}
tayltst_DEPENDENCIES=${filter %.a %.la,${tayltst_libs}}

