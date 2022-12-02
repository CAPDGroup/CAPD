include ${capdMake}/make/common_programs_makefile.mkf
ACLOCAL_AMFLAGS= -I ../../capdMake/ac_aux/m4

SUBDIRS=nonAutTest BasicCnTaylorTest diffIncl

if HAVE_BOOST
   SUBDIRS+=map dynset SolverGetStepTest ReversibleSolverTest ReversiblePoincareMapTest CnSolverTest MapMaskTest BasicSolverMaskTest BasicCurveMaskTest PoincareMapMaskTest
endif
