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
#m4_include([../../capdMake/make/configure_common_x.ac])


AC_CONFIG_FILES([Makefile])



AC_PROG_CXX
LT_INIT

capdkrak=yes
PKG_CHECK_MODULES([capdkrak], [capdkrak >= 1.0], [], [capdkrak=no])

AM_CONDITIONAL([HAVE_capdkrak], [test "x$capdkrak" = "xyes"])
if test "x$capdkrak" = "xyes"; then
  echo ${capdkrak_LIBS}
fi

AC_CONFIG_FILES([krakdemo/Makefile])
AC_CONFIG_FILES([configFileReader/Makefile])
AC_CONFIG_FILES([auxildemo/Makefile])
AC_CONFIG_FILES([logger/Makefile])

PKG_CHECK_MODULES([capdauxil], [capdauxil >= 1.0])

AC_CONFIG_HEADERS([config.h])
AC_OUTPUT
