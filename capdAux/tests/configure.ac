#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#AC_PREREQ([2.61])
m4_define([capdMake_PATH],[../../capdMake])
m4_include(capdMake_PATH[/make/capd_version.m4])
AC_INIT([capdAux-examples],[CAPD_VERSION_NUMBER],[http://capd.ii.uj.edu.pl],[],[http://capd.ii.uj.edu.pl])



AC_CONFIG_MACRO_DIRS([../../capdMake/ac_aux/m4 .]) # dot . for old autotools in Debian
AC_CONFIG_AUX_DIR(../../capdMake/ac_aux)
default_top_capd=${srcdir}/../..
top_capd_build=./../..
m4_include([../../capdMake/make/configure_common.ac])
AM_INIT_AUTOMAKE(capd_automake_init_options)

boost_optional=yes
m4_include([../../capdMake/make/configure_common_boost.ac])


AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([auxilTests/Makefile])


AC_PROG_CXX
LT_INIT

PKG_CHECK_MODULES([capdauxil], [capdauxil >= 1.0])
PKG_CHECK_MODULES([capdAux], [capdAux >= 1.0])


AM_COND_IF([HAVE_BOOST], [
  BOOST_TEST()
  BOOST_SYSTEM()
  BOOST_FILESYSTEM()
  BOOST_REGEX()
])

AC_CONFIG_HEADERS([config.h])
AC_OUTPUT
