#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#AC_PREREQ([2.61])
m4_define([capdMake_PATH],[../../capdMake])
m4_include(capdMake_PATH[/make/capd_version.m4])
AC_INIT([capdDynSys-examples],[CAPD_VERSION_NUMBER],[http://capd.ii.uj.edu.pl],[],[http://capd.ii.uj.edu.pl])

AC_CONFIG_MACRO_DIRS([../../capdMake/ac_aux/m4 .]) # dot . for old autotools in Debian
AC_CONFIG_AUX_DIR(../../capdMake/ac_aux)
default_top_capd=${srcdir}/../..
top_capd_build=./../..

m4_include([../../capdMake/make/configure_common.ac])
AM_INIT_AUTOMAKE(capd_automake_init_options)
m4_include([../../capdMake/make/configure_common_mpfr.ac])

PKG_CHECK_MODULES([capdAux], [capdAux >= 1.0])

HAVE_capdkrak=yes
PKG_CHECK_MODULES([capdkrak], [capdkrak >= 1.0], [], [HAVE_capdkrak=no])
AM_CONDITIONAL([HAVE_capdkrak], [test "x$HAVE_capdkrak" = "xyes"])

PKG_CHECK_MODULES([capdExt], [capdExt >= 1.0])
PKG_CHECK_MODULES([capdDynSys], [capdDynSys >= 1.0])
PKG_CHECK_MODULES([mpcapdDynSys], [mpcapdDynSys >= 1.0], [HAVE_mpcapdDynSys=yes], [AS_ECHO["No mpcapd, won't compile some examples"]])

AM_CONDITIONAL([HAVE_mpcapdDynSys], [test "x$HAVE_mpcapdDynSys" = "xyes"])
AM_CONDITIONAL([HAVE_OPENMP_CXX11], [false])


AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([maps/Makefile])
AC_CONFIG_FILES([odes/Makefile])
AC_CONFIG_FILES([odesrig/Makefile])
AC_CONFIG_FILES([cndemo/Makefile])
AC_CONFIG_FILES([coveringRelation/Makefile])
AC_CONFIG_FILES([diffIncldemo/Makefile])
AC_CONFIG_FILES([encloseTrajectory/Makefile])
AC_CONFIG_FILES([graphicalSet/Makefile])
AC_CONFIG_FILES([integrate/Makefile])
AC_CONFIG_FILES([multiPrecision/Makefile])
AC_CONFIG_FILES([newton/Makefile])
AC_CONFIG_FILES([poincare/Makefile])
AC_CONFIG_FILES([pointst/Makefile])
AC_CONFIG_FILES([taylor/Makefile])
AC_CONFIG_FILES([RosslerPeriodicOrbit/Makefile])
AC_CONFIG_FILES([RosslerChaoticDynamics/Makefile])
AC_CONFIG_FILES([LorenzPeriodicOrbit/Makefile])
AC_CONFIG_FILES([invariantSet/Makefile])



AC_PROG_CXX
LT_INIT

AC_CONFIG_HEADERS([config.h])
AC_OUTPUT
