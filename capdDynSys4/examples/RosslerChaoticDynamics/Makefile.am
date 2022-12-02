include ${capdMake}/make/common_programs_makefile.mkf

noinst_PROGRAMS=RosslerChaoticDynamics

RosslerChaoticDynamics_SOURCES=RosslerChaoticDynamics.cpp
RosslerChaoticDynamics_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS} ${CAPD_CXXFLAGS_CXX11}
RosslerChaoticDynamics_libs=${capdDynSys_LIBS}
RosslerChaoticDynamics_LDADD=${RosslerChaoticDynamics_libs}
RosslerChaoticDynamics_DEPENDENCIES=${filter %.a %.la,${RosslerChaoticDynamics_libs}}

EXTRA_DIST=attractor.png horseshoe.png RosslerChaoticDynamics.dox
