#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#AC_PREREQ([2.61])
m4_define([capdMake_PATH],[../capdMake])
m4_include(capdMake_PATH[/make/capd_version.m4])
AC_INIT([capdExt],[1.1],[http://capd.ii.uj.edu.pl],[],[http://capd.ii.uj.edu.pl])


AC_CONFIG_MACRO_DIRS([../capdMake/ac_aux/m4 .]) # dot . for old autotools in Debian
AC_CONFIG_AUX_DIR(../capdMake/ac_aux)
default_top_capd=${srcdir}/..
top_capd_build=./..
m4_include([../capdMake/m4/ax_mpi.m4])

m4_include([../capdMake/make/configure_common.ac])
AM_INIT_AUTOMAKE(capd_automake_init_options)

CAPD_CXX_FLAG_CHECK([-frounding-math], [frounding_math])
CAPD_CXX_FLAG_CHECK([-ffloat-store], [ffloat_store])

AC_ARG_WITH([cxsc],
    AS_HELP_STRING([--with-cxsc], [path to cxsc instalation. ]),
    [],
    [with_cxsc=no])

AC_ARG_WITH([filib],
    AS_HELP_STRING([--with-filib], [path to filib installation or internal to use filib version provided with CAPD. As a default we do not use filib.]),
    [],
    [with_filib=internal])

AS_IF([test "x$with_cxsc" = xyes],
    with_cxsc=/usr/local)

AS_IF([test "x$with_cxsc" != xno],
    [AC_CHECK_FILE([$with_cxsc/include/interval.hpp],[with_filib=no],[with_cxsc=no])])

AS_IF([test "x$with_cxsc" != xno],
      [AC_MSG_NOTICE([using CXSC intervals])
      AC_SUBST(filib_CFLAGS, ["-D__USE_CXSC__ -I$with_cxsc/include $CXX_FLAG_frounding_math"])
      AC_SUBST(filib_LIBS, ["-L$with_cxsc/lib -lcxsc"])
      AC_SUBST(filib_CFLAGS_INSTALL, ["-D__USE_CXSC__ -I$with_cxsc/include $CXX_FLAG_frounding_math"])
      AC_SUBST(filib_LIBS_INSTALL, ["-L$with_cxsc/lib -lcxsc"])])

AC_CANONICAL_HOST
echo HOST INFORMATION : $host
echo "$host_alias X $host_cpu X $host_vendor X  $host_os"

AS_IF([test "x$host_cpu" = xaarch64],
     [with_filib=no])

AS_IF([test "x$with_filib" = xcheck],
      [AC_CHECK_FILE([${srcdir}/../capdDynSys],[with_filib=internal],[with_filib=no])])

AS_IF([test "x$with_filib" = xyes],
       with_filib=internal)

AS_IF([test "x$with_filib" = xinternal],
      [AC_MSG_NOTICE([using internal filib])
      filib_CFLAGS="-D__USE_FILIB__ -I\${abs_top_srcdir}/filibsrc -I\${abs_top_builddir}/filibsrc $CXX_FLAG_frounding_math $CXX_FLAG_ffloat_store"
      filib_LIBS="\${abs_top_builddir}/filibsrc/libprim/libprim.la"
      AC_SUBST(filib_CFLAGS)
      AC_SUBST(filib_LIBS)
      AC_SUBST(filib_CFLAGS_INSTALL, ["-D__USE_FILIB__ -I\${prefix}/include $CXX_FLAG_frounding_math"])
      AC_SUBST(filib_LIBS_INSTALL, ["-L\${prefix}/lib -lprim"])])


AS_IF([test "x$with_filib" != xno && test "x$with_filib" != xinternal],
      [AC_MSG_NOTICE([using external filib])
      AC_SUBST(filib_CFLAGS, ["-D__USE_FILIB__ -I$with_filib/include $CXX_FLAG_frounding_math $CXX_FLAG_ffloat_store"])
      AC_SUBST(filib_LIBS, ["-L$with_filib/lib -lprim"])
      AC_SUBST(filib_CFLAGS_INSTALL, ["-D__USE_FILIB__ -I$with_filib/include $CXX_FLAG_frounding_math"])
      AC_SUBST(filib_LIBS_INSTALL, ["-L$with_filib/lib -lprim"])])

AS_IF([test "x$with_cxsc" = xno && test "x$with_filib" = xno],
      [AC_MSG_NOTICE([using CAPD intervals])
      AC_SUBST(filib_CFLAGS, [""])
      AC_SUBST(filib_LIBS, [""])
      AC_SUBST(filib_CFLAGS_INSTALL, [""])
      AC_SUBST(filib_LIBS_INSTALL, [""])])

AM_CONDITIONAL([WITH_FILIB_INTERNAL], [test "x$with_filib" = xinternal])



AC_SUBST(capdalglib_CFLAGS)
AC_CONFIG_FILES([${pkg_config_dev}/capdalglib.pc:src/capd/alglib/capdalglib-dev.pc.in])
AC_CONFIG_FILES([${pkg_config_install}/capdalglib.pc:src/capd/alglib/capdalglib-install.pc.in])
AC_SUBST(capdExt_CFLAGS)
AC_CONFIG_FILES([${pkg_config_dev}/capdExt.pc:src/capd/capdExt-dev.pc.in])
AC_CONFIG_FILES([${pkg_config_install}/capdExt.pc:src/capd/capdExt-install.pc.in])

AC_CONFIG_HEADERS([include/capd/config-capdExt.h])


AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([include/Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([src/capd/Makefile])
AC_CONFIG_FILES([src/capd/alglib/Makefile])

AS_IF([test "x$with_filib" = xinternal],
      [AC_CONFIG_SUBDIRS([filibsrc])])


#PKG_CHECK_MODULES([CAPD_PACKAGES], [capdAux >= 1.0])

# Checks for programs.
AC_PROG_CXX

# Checks for libraries.
LT_INIT

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
#AC_HEADER_STDBOOL

# Checks for library functions.

AC_OUTPUT
