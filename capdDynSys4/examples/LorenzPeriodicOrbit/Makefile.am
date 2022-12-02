include ${capdMake}/make/common_programs_makefile.mkf

noinst_PROGRAMS=LorenzPeriodicOrbit

LorenzPeriodicOrbit_SOURCES=LorenzPeriodicOrbit.cpp
LorenzPeriodicOrbit_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS} ${CAPD_CXXFLAGS_CXX11}
LorenzPeriodicOrbit_libs=${capdDynSys_LIBS}
LorenzPeriodicOrbit_LDADD=${LorenzPeriodicOrbit_libs}
LorenzPeriodicOrbit_DEPENDENCIES=${filter %.a %.la,${LorenzPeriodicOrbit_libs}}

EXTRA_DIST=LorenzPeriodicOrbits.dat LorenzPeriodicOrbit.dox
