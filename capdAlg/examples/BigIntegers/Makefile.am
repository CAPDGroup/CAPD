include ${capdMake}/make/common_programs_makefile.mkf


noinst_PROGRAMS=MpIntExample
MpIntExample_SOURCES=MpIntExample.cpp
MpIntExample_libs=${mpcapdAlg_LIBS} ${mpfr_LIBS} ${capdAux_LIBS}
MpIntExample_LDADD=${MpIntExample_libs}
MpIntExample_DEPENDENCIES=${filter %.a %.la,${MpIntExample_libs}}
MpIntExample_CXXFLAGS= ${AM_CXXFLAGS} ${mpcapdAlg_CFLAGS} ${mpfr_CFLAGS}  ${capdAux_CFLAGS}
