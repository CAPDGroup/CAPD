include ${capdMake}/make/common_programs_makefile.mkf

noinst_PROGRAMS=argdemo BuildInfo

argdemo_SOURCES=argdemo.cpp
argdemo_libs=${capdauxil_LIBS}
argdemo_LDADD=${argdemo_libs}
argdemo_DEPENDENCIES=${filter %.a %.la,${argdemo_libs}}
argdemo_CXXFLAGS= ${AM_CXXFLAGS} ${capdauxil_CFLAGS}

BuildInfo_SOURCES=BuildInfoMain.cpp
BuildInfo_libs=${capdauxil_LIBS}
BuildInfo_LDADD=${BuildInfo_libs}
BuildInfo_DEPENDENCIES=${filter %.a %.la,${BuildInfo_libs}}
BuildInfo_CXXFLAGS= ${AM_CXXFLAGS} ${capdauxil_CFLAGS}
