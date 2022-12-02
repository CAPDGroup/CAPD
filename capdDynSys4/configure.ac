#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#AC_PREREQ([2.61])
m4_define([capdMake_PATH],[../capdMake])
m4_include(capdMake_PATH[/make/capd_version.m4])
AC_INIT([capdDynSys],[CAPD_VERSION_NUMBER],[http://capd.ii.uj.edu.pl],[],[http://capd.ii.uj.edu.pl])

AC_CONFIG_MACRO_DIRS([../capdMake/ac_aux/m4 .]) # dot . for old autotools in Debian
AC_CONFIG_AUX_DIR(../capdMake/ac_aux)
default_top_capd=${srcdir}/..
top_capd_build=./..
m4_include([../capdMake/make/configure_common.ac])
AM_INIT_AUTOMAKE(capd_automake_init_options)
m4_include([../capdMake/make/configure_common_mpfr.ac])

# m4_include([../capdMake/make/configure_common_cxx11.ac])

m4_include([../capdMake/make/configure_common_cxx14.ac])

AC_SUBST(capdDynSys_CFLAGS)
AC_SUBST(capdDynSys_LIBS)
AC_CONFIG_FILES([${pkg_config_dev}/capdDynSys.pc:src/capd/capdDynSys-dev.pc.in])
AC_CONFIG_FILES([${pkg_config_install}/capdDynSys.pc:src/capd/capdDynSys-install.pc.in])
AC_SUBST(capdautodiff_CFLAGS)
AC_SUBST(capdautodiff_LIBS)
AC_CONFIG_FILES([${pkg_config_dev}/capdautodiff.pc:src/capd/autodiff/capdautodiff-dev.pc.in])
AC_CONFIG_FILES([${pkg_config_install}/capdautodiff.pc:src/capd/autodiff/capdautodiff-install.pc.in])
AC_SUBST(capdcovrel_CFLAGS)
AC_SUBST(capdcovrel_LIBS)
AC_CONFIG_FILES([${pkg_config_dev}/capdcovrel.pc:src/capd/covrel/capdcovrel-dev.pc.in])
AC_CONFIG_FILES([${pkg_config_install}/capdcovrel.pc:src/capd/covrel/capdcovrel-install.pc.in])
AC_SUBST(capddiffIncl_CFLAGS)
AC_SUBST(capddiffIncl_LIBS)
AC_CONFIG_FILES([${pkg_config_dev}/capddiffIncl.pc:src/capd/diffIncl/capddiffIncl-dev.pc.in])
AC_CONFIG_FILES([${pkg_config_install}/capddiffIncl.pc:src/capd/diffIncl/capddiffIncl-install.pc.in])
AC_SUBST(capddynsysm_CFLAGS)
AC_SUBST(capddynsysm_LIBS)
AC_CONFIG_FILES([${pkg_config_dev}/capddynsysm.pc:src/capd/dynsys/capddynsysm-dev.pc.in])
AC_CONFIG_FILES([${pkg_config_install}/capddynsysm.pc:src/capd/dynsys/capddynsysm-install.pc.in])
AC_SUBST(capdgeomset_CFLAGS)
AC_SUBST(capdgeomset_LIBS)
AC_CONFIG_FILES([${pkg_config_dev}/capdgeomset.pc:src/capd/geomset/capdgeomset-dev.pc.in])
AC_CONFIG_FILES([${pkg_config_install}/capdgeomset.pc:src/capd/geomset/capdgeomset-install.pc.in])
AC_SUBST(capdmap_CFLAGS)
AC_SUBST(capdmap_LIBS)
AC_CONFIG_FILES([${pkg_config_dev}/capdmap.pc:src/capd/map/capdmap-dev.pc.in])
AC_CONFIG_FILES([${pkg_config_install}/capdmap.pc:src/capd/map/capdmap-install.pc.in])
AC_SUBST(capddiffAlgebra_CFLAGS)
AC_SUBST(capddiffAlgebra_LIBS)
AC_CONFIG_FILES([${pkg_config_dev}/capddiffAlgebra.pc:src/capd/diffAlgebra/capddiffAlgebra-dev.pc.in])
AC_CONFIG_FILES([${pkg_config_install}/capddiffAlgebra.pc:src/capd/diffAlgebra/capddiffAlgebra-install.pc.in])
AC_SUBST(capddynset_CFLAGS)
AC_SUBST(capddynset_LIBS)
AC_CONFIG_FILES([${pkg_config_dev}/capddynset.pc:src/capd/dynset/capddynset-dev.pc.in])
AC_CONFIG_FILES([${pkg_config_install}/capddynset.pc:src/capd/dynset/capddynset-install.pc.in])
AC_SUBST(capdpoincare_CFLAGS)
AC_SUBST(capdpoincare_LIBS)
AC_CONFIG_FILES([${pkg_config_dev}/capdpoincare.pc:src/capd/poincare/capdpoincare-dev.pc.in])
AC_CONFIG_FILES([${pkg_config_install}/capdpoincare.pc:src/capd/poincare/capdpoincare-install.pc.in])
AC_SUBST(capdpdes_CFLAGS)
AC_SUBST(capdpdes_LIBS)
AC_CONFIG_FILES([${pkg_config_dev}/capdpdes.pc:src/capd/pdes/capdpdes-dev.pc.in])
AC_CONFIG_FILES([${pkg_config_install}/capdpdes.pc:src/capd/pdes/capdpdes-install.pc.in])



AC_CONFIG_HEADERS([include/capd/config-capdDynSys.h])


AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([include/Makefile])
AC_CONFIG_FILES([src/Makefile])

AC_CONFIG_FILES([src/capd/Makefile])
AC_CONFIG_FILES([src/capd/autodiff/Makefile])
AC_CONFIG_FILES([src/capd/covrel/Makefile])
AC_CONFIG_FILES([src/capd/diffIncl/Makefile])
AC_CONFIG_FILES([src/capd/dynsys/Makefile])
AC_CONFIG_FILES([src/capd/geomset/Makefile])
AC_CONFIG_FILES([src/capd/map/Makefile])
AC_CONFIG_FILES([src/capd/diffAlgebra/Makefile])
AC_CONFIG_FILES([src/capd/dynset/Makefile])
AC_CONFIG_FILES([src/capd/poincare/Makefile])
AC_CONFIG_FILES([src/capd/pdes/Makefile])

PKG_CHECK_EXISTS([mpcapdAlg], [HAVE_MPCAPDALG=yes])
AM_CONDITIONAL([HAVE_MPCAPDALG],[test "x$HAVE_MPCAPDALG" = "xyes"])

AS_IF([test "x$HAVE_MPCAPDALG" = "xyes"],
   [PKG_CHECK_MODULES([mpcapdAlg], [mpcapdAlg >= 1.0])
   AC_SUBST(mpcapdDynSys_CFLAGS)
   AC_SUBST(mpcapdDynSys_LIBS)
    AC_CONFIG_FILES([${pkg_config_dev}/mpcapdDynSys.pc:src/mpcapd/mpcapdDynSys-dev.pc.in])
    AC_CONFIG_FILES([${pkg_config_install}/mpcapdDynSys.pc:src/mpcapd/mpcapdDynSys-install.pc.in])

    AC_SUBST(mpcapdautodiff_CFLAGS, ["$mpfr_CFLAGS -frounding-math"])
    AC_SUBST(mpcapdautodiff_LIBS, ["$mpfr_LIBS"])
    AC_CONFIG_FILES([${pkg_config_dev}/mpcapdautodiff.pc:src/mpcapd/autodiff/mpcapdautodiff-dev.pc.in])
    AC_CONFIG_FILES([${pkg_config_install}/mpcapdautodiff.pc:src/mpcapd/autodiff/mpcapdautodiff-install.pc.in])


    AC_SUBST(mpcapddiffIncl_CFLAGS, ["$mpfr_CFLAGS -frounding-math"])
    AC_SUBST(mpcapddiffIncl_LIBS, ["$mpfr_LIBS"])
    AC_CONFIG_FILES([${pkg_config_dev}/mpcapddiffIncl.pc:src/mpcapd/diffIncl/mpcapddiffIncl-dev.pc.in])
    AC_CONFIG_FILES([${pkg_config_install}/mpcapddiffIncl.pc:src/mpcapd/diffIncl/mpcapddiffIncl-install.pc.in])

    AC_SUBST(mpcapddynsysm_CFLAGS, ["$mpfr_CFLAGS -frounding-math"])
    AC_SUBST(mpcapddynsysm_LIBS, ["$mpfr_LIBS"])
    AC_CONFIG_FILES([${pkg_config_dev}/mpcapddynsysm.pc:src/mpcapd/dynsys/mpcapddynsysm-dev.pc.in])
    AC_CONFIG_FILES([${pkg_config_install}/mpcapddynsysm.pc:src/mpcapd/dynsys/mpcapddynsysm-install.pc.in])

    AC_SUBST(mpcapdgeomset_CFLAGS, ["$mpfr_CFLAGS -frounding-math"])
    AC_SUBST(mpcapdgeomset_LIBS, ["$mpfr_LIBS"])
    AC_CONFIG_FILES([${pkg_config_dev}/mpcapdgeomset.pc:src/mpcapd/geomset/mpcapdgeomset-dev.pc.in])
    AC_CONFIG_FILES([${pkg_config_install}/mpcapdgeomset.pc:src/mpcapd/geomset/mpcapdgeomset-install.pc.in])

    AC_SUBST(mpcapdmap_CFLAGS, ["$mpfr_CFLAGS -frounding-math"])
    AC_SUBST(mpcapdmap_LIBS, ["$mpfr_LIBS"])
    AC_CONFIG_FILES([${pkg_config_dev}/mpcapdmap.pc:src/mpcapd/map/mpcapdmap-dev.pc.in])
    AC_CONFIG_FILES([${pkg_config_install}/mpcapdmap.pc:src/mpcapd/map/mpcapdmap-install.pc.in])

    AC_SUBST(mpcapddiffAlgebra_CFLAGS, ["$mpfr_CFLAGS -frounding-math"])
    AC_SUBST(mpcapddiffAlgebra_LIBS, ["$mpfr_LIBS"])
    AC_CONFIG_FILES([${pkg_config_dev}/mpcapddiffAlgebra.pc:src/mpcapd/diffAlgebra/mpcapddiffAlgebra-dev.pc.in])
    AC_CONFIG_FILES([${pkg_config_install}/mpcapddiffAlgebra.pc:src/mpcapd/diffAlgebra/mpcapddiffAlgebra-install.pc.in])

    AC_SUBST(mpcapddynset_CFLAGS, ["$mpfr_CFLAGS -frounding-math"])
    AC_SUBST(mpcapddynset_LIBS, ["$mpfr_LIBS"])
    AC_CONFIG_FILES([${pkg_config_dev}/mpcapddynset.pc:src/mpcapd/dynset/mpcapddynset-dev.pc.in])
    AC_CONFIG_FILES([${pkg_config_install}/mpcapddynset.pc:src/mpcapd/dynset/mpcapddynset-install.pc.in])

    AC_SUBST(mpcapdpoincare_CFLAGS, ["$mpfr_CFLAGS -frounding-math"])
    AC_SUBST(mpcapdpoincare_LIBS, ["$mpfr_LIBS"])
    AC_CONFIG_FILES([${pkg_config_dev}/mpcapdpoincare.pc:src/mpcapd/poincare/mpcapdpoincare-dev.pc.in])
    AC_CONFIG_FILES([${pkg_config_install}/mpcapdpoincare.pc:src/mpcapd/poincare/mpcapdpoincare-install.pc.in])

    AC_CONFIG_FILES([src/mpcapd/Makefile])
    AC_CONFIG_FILES([src/mpcapd/autodiff/Makefile])
    AC_CONFIG_FILES([src/mpcapd/diffAlgebra/Makefile])
    AC_CONFIG_FILES([src/mpcapd/diffIncl/Makefile])
    AC_CONFIG_FILES([src/mpcapd/dynset/Makefile])
    AC_CONFIG_FILES([src/mpcapd/dynsys/Makefile])
    AC_CONFIG_FILES([src/mpcapd/geomset/Makefile])
    AC_CONFIG_FILES([src/mpcapd/map/Makefile])
    AC_CONFIG_FILES([src/mpcapd/poincare/Makefile])],
   [rm -f ${top_capd_build}/lib/pkgconfig-dev/mpcapdDynSys.pc
    rm -f ${top_capd_build}/lib/pkgconfig-dev/mpcapdautodiff.pc
    rm -f ${top_capd_build}/lib/pkgconfig-dev/mpcapddiffIncl.pc
    rm -f ${top_capd_build}/lib/pkgconfig-dev/mpcapddynsysm.pc
    rm -f ${top_capd_build}/lib/pkgconfig-dev/mpcapdgeomset.pc
    rm -f ${top_capd_build}/lib/pkgconfig-dev/mpcapdmap.pc
    rm -f ${top_capd_build}/lib/pkgconfig-dev/mpcapddiffAlgebra.pc
    rm -f ${top_capd_build}/lib/pkgconfig-dev/mpcapdpoincare.pc
    rm -f ${top_capd_build}/lib/pkgconfig-dev/mpcapddynset.pc])


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
PKG_CHECK_MODULES([capdAlg], [capdAlg >= 1.0])
PKG_CHECK_MODULES([capdExt], [capdExt >= 1.0])


AC_OUTPUT
