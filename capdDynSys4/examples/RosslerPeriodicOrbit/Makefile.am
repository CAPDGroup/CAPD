include ${capdMake}/make/common_programs_makefile.mkf

noinst_PROGRAMS=RosslerPeriodicOrbit

RosslerPeriodicOrbit_SOURCES=RosslerPeriodicOrbit.cpp
RosslerPeriodicOrbit_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS} ${CAPD_CXXFLAGS_CXX11}
RosslerPeriodicOrbit_libs=${capdDynSys_LIBS}
RosslerPeriodicOrbit_LDADD=${RosslerPeriodicOrbit_libs}
RosslerPeriodicOrbit_DEPENDENCIES=${filter %.a %.la,${RosslerPeriodicOrbit_libs}}

EXTRA_DIST=Rossler-3per-orbits.png RosslerPeriodicOrbit.dox
