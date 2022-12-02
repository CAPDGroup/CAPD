include ${capdMake}/make/common_programs_makefile.mkf

noinst_PROGRAMS=logger
logger_SOURCES=logger.cpp Example.h
logger_libs=${capdauxil_LIBS}
logger_LDADD=${logger_libs}
logger_DEPENDENCIES=${filter %.a %.la,${logger_libs}}
logger_CXXFLAGS= ${AM_CXXFLAGS} ${capdauxil_CFLAGS}
