#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#AC_PREREQ([2.61])
m4_define([capdMake_PATH],[../capdMake])
m4_include(capdMake_PATH[/make/capd_version.m4])
AC_INIT([capdAlg],[CAPD_VERSION_NUMBER],[http://capd.ii.uj.edu.pl],[],[http://capd.ii.uj.edu.pl])

AC_CONFIG_MACRO_DIRS([../capdMake/ac_aux/m4 .]) # dot . for old autotools in Debian
AC_CONFIG_AUX_DIR(../capdMake/ac_aux)
default_top_capd=${srcdir}/..
top_capd_build=./..
m4_include([../capdMake/make/configure_common.ac])
AM_INIT_AUTOMAKE(capd_automake_init_options)
m4_include([../capdMake/make/configure_common_mpfr.ac])
m4_include([../capdMake/make/configure_common_pari.ac])
#m4_include([../capdMake/make/configure_common_cxx11.ac])

CAPD_AM_CXXFLAGS="${CAPD_CXXFLAGS_CXX11} ${CAPD_AM_CXXFLAGS}" # only for capd compilation, do not have to add it to capd<module>_CFLAGS.
AC_SUBST([CAPD_AM_CXXFLAGS])

AS_IF([test "x$HAVE_PARI" = "xyes"],
[
  capdmatrixAlgorithms_LIBS="$capdmatrixAlgorithms_LIBS -lpari"
])

mpcapdmatrixAlgorithms_LIBS="$capdmatrixAlgorithms_LIBS"

CAPD_CXX_FLAG_CHECK([-frounding-math], [frounding_math])

AS_IF([test "x$HAVE_LIBmpfr" = "xyes"],
	    [
	    AC_DEFINE([HAVE_MPCAPD_ALG], [1], [Have mpcapdAlg])
            mpcapdAlg_CFLAGS="${mpcapdAlg_CFLAGS} -DHAVE_MPCAPD_ALG"

	    AC_SUBST(mpcapdAlg_CFLAGS)
	    AC_SUBST(mpcapdAlg_LIBS)
	    AC_CONFIG_FILES([${pkg_config_dev}/mpcapdAlg.pc:src/mpcapd/mpcapdAlg-dev.pc.in])
	    AC_CONFIG_FILES([${pkg_config_install}/mpcapdAlg.pc:src/mpcapd/mpcapdAlg-install.pc.in])
	    AC_SUBST(mpcapdintervals_CFLAGS, ["$mpfr_CFLAGS $CXX_FLAG_frounding_math"])
	    AC_SUBST(mpcapdintervals_LIBS, ["$mpfr_LIBS"])
	    AC_CONFIG_FILES([${pkg_config_dev}/mpcapdintervals.pc:src/mpcapd/intervals/mpcapdintervals-dev.pc.in])
	    AC_CONFIG_FILES([${pkg_config_install}/mpcapdintervals.pc:src/mpcapd/intervals/mpcapdintervals-install.pc.in])
	    AC_SUBST(mpcapdmatrixAlgorithms_CFLAGS, ["$mpfr_CFLAGS $CXX_FLAG_frounding_math"])
	    AC_SUBST(mpcapdmatrixAlgorithms_LIBS, ["$mpfr_LIBS"])
	    AC_CONFIG_FILES([${pkg_config_dev}/mpcapdmatrixAlgorithms.pc:src/mpcapd/matrixAlgorithms/mpcapdmatrixAlgorithms-dev.pc.in])
	    AC_CONFIG_FILES([${pkg_config_install}/mpcapdmatrixAlgorithms.pc:src/mpcapd/matrixAlgorithms/mpcapdmatrixAlgorithms-install.pc.in])
	    AC_SUBST(mpcapdvectalg_CFLAGS, ["$mpfr_CFLAGS $CXX_FLAG_frounding_math"])
	    AC_SUBST(mpcapdvectalg_LIBS, ["$mpfr_LIBS"])
	    AC_CONFIG_FILES([${pkg_config_dev}/mpcapdvectalg.pc:src/mpcapd/vectalg/mpcapdvectalg-dev.pc.in])
	    AC_CONFIG_FILES([${pkg_config_install}/mpcapdvectalg.pc:src/mpcapd/vectalg/mpcapdvectalg-install.pc.in])
	    echo
	    echo Multiprecision support enabled
	    echo
	    ],
	    [
	    rm -f ${top_capd_build}/lib/pkgconfig-dev/mpcapdAlg.pc
	    rm -f ${top_capd_build}/lib/pkgconfig-dev/mpcapdintervals.pc
	    rm -f ${top_capd_build}/lib/pkgconfig-dev/mpcapdmatrixAlgorithms.pc
	    rm -f ${top_capd_build}/lib/pkgconfig-dev/mpcapdvectalg.pc
	    ]
)


AC_SUBST(capdAlg_CFLAGS)
AC_SUBST(capdAlg_LIBS)
AC_CONFIG_FILES([${pkg_config_dev}/capdAlg.pc:src/capd/capdAlg-dev.pc.in])
AC_CONFIG_FILES([${pkg_config_install}/capdAlg.pc:src/capd/capdAlg-install.pc.in])
AC_SUBST(capdbasicalg_CFLAGS, ["$CXX_FLAG_frounding_math"])
AC_SUBST(capdbasicalg_LIBS)
AC_CONFIG_FILES([${pkg_config_dev}/capdbasicalg.pc:src/capd/basicalg/capdbasicalg-dev.pc.in])
AC_CONFIG_FILES([${pkg_config_install}/capdbasicalg.pc:src/capd/basicalg/capdbasicalg-install.pc.in])
AC_SUBST(capdintervals_CFLAGS, ["$CXX_FLAG_frounding_math"])
AC_SUBST(capdintervals_LIBS)
AC_CONFIG_FILES([${pkg_config_dev}/capdintervals.pc:src/capd/intervals/capdintervals-dev.pc.in])
AC_CONFIG_FILES([${pkg_config_install}/capdintervals.pc:src/capd/intervals/capdintervals-install.pc.in])
AC_SUBST(capdmatrixAlgorithms_CFLAGS, ["$CXX_FLAG_frounding_math"])
AC_SUBST(capdmatrixAlgorithms_LIBS)
AC_CONFIG_FILES([${pkg_config_dev}/capdmatrixAlgorithms.pc:src/capd/matrixAlgorithms/capdmatrixAlgorithms-dev.pc.in])
AC_CONFIG_FILES([${pkg_config_install}/capdmatrixAlgorithms.pc:src/capd/matrixAlgorithms/capdmatrixAlgorithms-install.pc.in])

AC_SUBST(capdrings_CFLAGS)
AC_SUBST(capdrings_LIBS)
AC_CONFIG_FILES([${pkg_config_dev}/capdrings.pc:src/capd/rings/capdrings-dev.pc.in])
AC_CONFIG_FILES([${pkg_config_install}/capdrings.pc:src/capd/rings/capdrings-install.pc.in])

AC_SUBST(capdrounding_CFLAGS, ["$CXX_FLAG_frounding_math"])
AC_SUBST(capdrounding_LIBS)
AC_CONFIG_FILES([${pkg_config_dev}/capdrounding.pc:src/capd/rounding/capdrounding-dev.pc.in])
AC_CONFIG_FILES([${pkg_config_install}/capdrounding.pc:src/capd/rounding/capdrounding-install.pc.in])

AC_SUBST(capdvectalg_CFLAGS, ["$CXX_FLAG_frounding_math"])
AC_SUBST(capdvectalg_LIBS)
AC_CONFIG_FILES([${pkg_config_dev}/capdvectalg.pc:src/capd/vectalg/capdvectalg-dev.pc.in])
AC_CONFIG_FILES([${pkg_config_install}/capdvectalg.pc:src/capd/vectalg/capdvectalg-install.pc.in])

AC_CONFIG_HEADERS([include/capd/config-capdAlg.h])


AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([include/Makefile])
AC_CONFIG_FILES([src/Makefile])

AC_CONFIG_FILES([src/capd/Makefile])
AC_CONFIG_FILES([src/capd/basicalg/Makefile])
AC_CONFIG_FILES([src/capd/intervals/Makefile])
AC_CONFIG_FILES([src/capd/matrixAlgorithms/Makefile])
AC_CONFIG_FILES([src/capd/rings/Makefile])
AC_CONFIG_FILES([src/capd/rounding/Makefile])
AC_CONFIG_FILES([src/capd/vectalg/Makefile])

AC_CONFIG_FILES([src/mpcapd/Makefile])
AC_CONFIG_FILES([src/mpcapd/intervals/Makefile])
AC_CONFIG_FILES([src/mpcapd/matrixAlgorithms/Makefile])
AC_CONFIG_FILES([src/mpcapd/vectalg/Makefile])


AC_CONFIG_SUBDIRS([examples])
AC_CONFIG_SUBDIRS([tests])


# Checks for programs.
AC_PROG_CXX

# Checks for libraries.
LT_INIT

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
#AC_HEADER_STDBOOL

# Checks for library functions.
PKG_CHECK_MODULES([capdAux], [capdAux >= 1.0])
PKG_CHECK_MODULES([capdExt], [capdExt >= 1.0])

AC_OUTPUT
