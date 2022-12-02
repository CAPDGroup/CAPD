include ${capdMake}/make/common_programs_makefile.mkf
ACLOCAL_AMFLAGS= -I ../../capdMake/ac_aux/m4

SUBDIRS=intervalTest  intvtst  matrixAlgorithms  vecttst basicalg complex

if HAVE_LIBmpfr
   SUBDIRS += mptests
endif
