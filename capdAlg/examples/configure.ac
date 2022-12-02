#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#AC_PREREQ([2.61])
m4_define([capdMake_PATH],[../../capdMake])
m4_include(capdMake_PATH[/make/capd_version.m4])
AC_INIT([capdAlg-examples],[CAPD_VERSION_NUMBER],[http://capd.ii.uj.edu.pl],[],[http://capd.ii.uj.edu.pl])

AC_CONFIG_MACRO_DIRS([../../capdMake/ac_aux/m4 .]) # dot . for old autotools in Debian
AC_CONFIG_AUX_DIR(../../capdMake/ac_aux)
default_top_capd=${srcdir}/../..
top_capd_build=./../..

m4_include([../../capdMake/make/configure_common.ac])
AM_INIT_AUTOMAKE(capd_automake_init_options)
m4_include([../../capdMake/make/configure_common_mpfr.ac])

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([alglibdemo/Makefile])
#AC_CONFIG_FILES([BigIntegers/Makefile])
AC_CONFIG_FILES([intervals/Makefile])
AC_CONFIG_FILES([matrixAlgorithms/Makefile])
AC_CONFIG_FILES([vectalg/Makefile])
AC_CONFIG_FILES([TexWriter/Makefile])

AC_PROG_CXX
LT_INIT

PKG_CHECK_MODULES([capdbasicalg], [capdbasicalg >= 1.0])
PKG_CHECK_MODULES([capdintervals], [capdintervals >= 1.0])
PKG_CHECK_MODULES([capdmatrixAlgorithms], [capdmatrixAlgorithms >= 1.0])
PKG_CHECK_MODULES([capdrounding], [capdrounding >= 1.0])
PKG_CHECK_MODULES([capdvectalg], [capdvectalg >= 1.0])

PKG_CHECK_MODULES([capdAux], [capdAux >= 1.0])
PKG_CHECK_MODULES([capdExt], [capdExt >= 1.0])

PKG_CHECK_MODULES([mpcapdAlg], [mpcapdAlg >= 1.0], [], [AC_MSG_WARN("No mpcapd")])

AC_CONFIG_HEADERS([config.h])
AC_OUTPUT
