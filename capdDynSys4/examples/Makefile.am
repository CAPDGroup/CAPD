include ${capdMake}/make/common_programs_makefile.mkf
ACLOCAL_AMFLAGS= -I ../../capdMake/ac_aux/m4

SUBDIRS=cndemo diffIncldemo encloseTrajectory integrate newton coveringRelation maps odes odesrig poincare
EXTRA_DIST=examples.dox

if HAVE_capdkrak
   SUBDIRS+= graphicalSet taylor pointst
endif

if HAVE_mpcapdDynSys
   SUBDIRS += multiPrecision
endif


   SUBDIRS += RosslerPeriodicOrbit LorenzPeriodicOrbit RosslerChaoticDynamics

# the flag is not set in configure.ac
if HAVE_OPENMP_CXX11
   SUBDIRS += invariantSet
endif
