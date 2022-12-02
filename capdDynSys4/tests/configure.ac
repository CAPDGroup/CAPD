#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#AC_PREREQ([2.61])
m4_define([capdMake_PATH],[../../capdMake])
m4_include(capdMake_PATH[/make/capd_version.m4])
AC_INIT([capdDynSys-tests],[CAPD_VERSION_NUMBER],[wilczak@ii.uj.edu.pl])



AC_CONFIG_MACRO_DIRS([../../capdMake/ac_aux/m4 .]) # dot . for old autotools in Debian
AC_CONFIG_AUX_DIR(../../capdMake/ac_aux)
default_top_capd=${srcdir}/../..
top_capd_build=./../..
m4_include([../../capdMake/make/configure_common.ac])
AM_INIT_AUTOMAKE(capd_automake_init_options)
boost_optional=yes
m4_include([../../capdMake/make/configure_common_boost.ac])

AC_CONFIG_FILES([BasicCnTaylorTest/Makefile])
AC_CONFIG_FILES([BasicCurveMaskTest/Makefile])
AC_CONFIG_FILES([BasicSolverMaskTest/Makefile])
AC_CONFIG_FILES([CnSolverTest/Makefile])
AC_CONFIG_FILES([diffIncl/Makefile])
AC_CONFIG_FILES([dynset/Makefile])
AC_CONFIG_FILES([map/Makefile])
AC_CONFIG_FILES([MapMaskTest/Makefile])
AC_CONFIG_FILES([nonAutTest/Makefile])
AC_CONFIG_FILES([PoincareMapMaskTest/Makefile])
AC_CONFIG_FILES([ReversiblePoincareMapTest/Makefile])
AC_CONFIG_FILES([ReversibleSolverTest/Makefile])
AC_CONFIG_FILES([SolverGetStepTest/Makefile])
AC_CONFIG_FILES([Makefile])

AC_PROG_CXX
LT_INIT

PKG_CHECK_MODULES([capdAux], [capdAux >= 1.0])
PKG_CHECK_MODULES([capdExt], [capdExt >= 1.0])
PKG_CHECK_MODULES([capdDynSys], [capdDynSys >= 1.0])


HAVE_capdkrak=yes
PKG_CHECK_MODULES([capdkrak], [capdkrak >= 1.0], [], [HAVE_capdkrak=no])
AM_CONDITIONAL([HAVE_capdkrak], [test "x$HAVE_capdkrak" = "xyes"])


AM_COND_IF([HAVE_BOOST],[
    BOOST_TEST()
])

AC_CONFIG_HEADERS([config.h])
AC_OUTPUT
