#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.


#AC_PREREQ([2.61])
m4_define([capdMake_PATH],[../capdMake])
m4_include(capdMake_PATH[/make/capd_version.m4])
AC_INIT([capdAux],[CAPD_VERSION_NUMBER],[http://capd.ii.uj.edu.pl],[],[http://capd.ii.uj.edu.pl])

AC_CONFIG_MACRO_DIRS([../capdMake/ac_aux/m4 .]) # dot . for old autotools in Debian
AC_CONFIG_AUX_DIR(../capdMake/ac_aux)
default_top_capd=${srcdir}/..
top_capd_build=./..


m4_include([../capdMake/make/configure_common.ac])
AM_INIT_AUTOMAKE(capd_automake_init_options)
m4_include([../capdMake/make/configure_common_x.ac])

boost_optional=yes
m4_include([../capdMake/make/configure_common_boost.ac])


AC_CANONICAL_HOST
echo HOST INFORMATION : $host
echo "$host_alias X $host_cpu X $host_vendor X  $host_os"

AC_PROG_CXX
LT_INIT


AM_COND_IF([HAVE_BOOST], [
  BOOST_SYSTEM()
  BOOST_FILESYSTEM()
  BOOST_REGEX()
])

AC_SUBST([DEF_SVN_REVISION_1], [$SVN_REVISION_1])
AC_DEFINE_UNQUOTED([DEF_SVN_REVISION_1], ["$SVN_REVISION_1"], [DEF_SVN_REVISION_1])
AC_SUBST([DEF_SVN_REVISION_2], [$SVN_REVISION_2])
AC_DEFINE_UNQUOTED([DEF_SVN_REVISION_2], ["$SVN_REVISION_2"], [DEF_SVN_REVISION_2])
AC_SUBST([DEF_SVN_REVISION_3], [$SVN_REVISION_3])
AC_DEFINE_UNQUOTED([DEF_SVN_REVISION_3], ["$SVN_REVISION_3"], [DEF_SVN_REVISION_3])
AC_SUBST([DEF_SVN_REVISION_4], [$SVN_REVISION_4])
AC_DEFINE_UNQUOTED([DEF_SVN_REVISION_4], ["$SVN_REVISION_4"], [DEF_SVN_REVISION_4])
AC_SUBST([DEF_SVN_URL_1], [$SVN_URL_1])
AC_DEFINE_UNQUOTED([DEF_SVN_URL_1], ["$SVN_URL_1"], [DEF_SVN_URL_1])
AC_SUBST([DEF_SVN_URL_2], [$SVN_URL_2])
AC_DEFINE_UNQUOTED([DEF_SVN_URL_2], ["$SVN_URL_2"], [DEF_SVN_URL_2])
AC_SUBST([DEF_SVN_URL_3], [$SVN_URL_3])
AC_DEFINE_UNQUOTED([DEF_SVN_URL_3], ["$SVN_URL_3"], [DEF_SVN_URL_3])
AC_SUBST([DEF_SVN_URL_4], [$SVN_URL_4])
AC_DEFINE_UNQUOTED([DEF_SVN_URL_4], ["$SVN_URL_4"], [DEF_SVN_URL_4])
AC_SUBST([DEF_BUILD_TAG], [$BUILD_TAG])
AC_DEFINE_UNQUOTED([DEF_BUILD_TAG], ["$BUILD_TAG"], [DEF_BUILD_TAG])
AC_SUBST([DEF_JOB_NAME], [$JOB_NAME])
AC_DEFINE_UNQUOTED([DEF_JOB_NAME], ["$JOB_NAME"], [DEF_JOB_NAME])
AC_SUBST([DEF_BUILD_DISPLAY_NAME], [$BUILD_DISPLAY_NAME])
AC_DEFINE_UNQUOTED([DEF_BUILD_DISPLAY_NAME], ["$BUILD_DISPLAY_NAME"], [DEF_BUILD_DISPLAY_NAME])
AC_SUBST([DEF_BUILD_ID], [$BUILD_ID])
AC_DEFINE_UNQUOTED([DEF_BUILD_ID], ["$BUILD_ID"], [DEF_BUILD_ID])
AC_SUBST([DEF_NODE_NAME], [$NODE_NAME])
AC_DEFINE_UNQUOTED([DEF_NODE_NAME], ["$NODE_NAME"], [DEF_NODE_NAME])
AC_SUBST([DEF_BUILD_NUMBER], [$BUILD_NUMBER])
AC_DEFINE_UNQUOTED([DEF_BUILD_NUMBER], ["$BUILD_NUMBER"], [DEF_BUILD_NUMBER])
AC_SUBST([DEF_VERSION], [$VERSION])
AC_DEFINE_UNQUOTED([DEF_VERSION], ["$VERSION"], [DEF_VERSION])
AC_SUBST([DEF_BUILD_DATE], [$(date +%Y%m%d-%H%M)])
AC_DEFINE_UNQUOTED([DEF_BUILD_DATE], ["$(date +%Y%m%d-%H%M)"], [DEF_BUILD_DATE])
AC_SUBST([DEF_LIBRARY_VERSION], ["$CAPD_VERSION_INFO"])
AC_DEFINE_UNQUOTED([DEF_LIBRARY_VERSION], ["$CAPD_VERSION_INFO"], [DEF_LIBRARY_VERSION])


PKG_CHECK_MODULES([liblog4cxx], [liblog4cxx >= 0.10.0], [AC_DEFINE([HAVE_LOG4CXX], [1], [Use log4cxx])],
				 [
				 AC_MSG_WARN([Cannot find log4cxx])
				 ])

CAPD_MODULES="capdauxil "

capdauxil_CFLAGS="${capdauxil_CFLAGS} ${liblog4cxx_CFLAGS}"
capdauxil_LIBS="${capdauxil_LIBS} ${liblog4cxx_LIBS}"

if test x"$HAVE_BOOST_system" = xyes && test x"$HAVE_BOOST_filesystem" = xyes && test x"$HAVE_BOOST_regex" = xyes; then
   capdauxil_CFLAGS="${capdauxil_CFLAGS} ${BOOST_CPPFLAGS}"
   capdauxil_LIBS="${capdauxil_LIBS} ${BOOST_FILESYSTEM_LIBS} ${BOOST_SYSTEM_LIBS} ${BOOST_REGEX_LIBS} ${BOOST_FILESYSTEM_LDFLAGS} ${BOOST_SYSTEM_LDFLAGS} ${BOOST_REGEX_LDFLAGS}"
else
  AC_MSG_WARN([Cannot find Boost: system or filesystem or regex])
fi


AC_SUBST(capdauxil_CFLAGS)
AC_SUBST(capdauxil_LIBS)
AC_CONFIG_FILES([${pkg_config_dev}/capdauxil.pc:src/capd/auxil/capdauxil-dev.pc.in])
AC_CONFIG_FILES([${pkg_config_install}/capdauxil.pc:src/capd/auxil/capdauxil-install.pc.in])

rm -f ${top_capd_build}/lib/pkgconfig-dev/capdkrak.pc


capdkrak_CFLAGS=""
capdkrak_LIBS=""
capdkrak_REQUIRES=""

AS_IF([test "x$HAVE_GUI" = "xyes"],
   [
#capdkrak isn't default, returned only by capd-config-gui or manual
#   CAPD_MODULES="$CAPD_MODULES capdkrak"

    capdkrak_CFLAGS=${GUI_FLAGS}
    capdkrak_LIBS=${GUI_LIBS}
    capdkrak_REQUIRES=${GUI_REQUIRES}

    AC_SUBST(capdkrak_CFLAGS)
    AC_SUBST(capdkrak_LIBS)
    AC_SUBST(capdkrak_REQUIRES)

    AC_CONFIG_FILES([${pkg_config_dev}/capdkrak.pc:src/capd/krak/capdkrak-dev.pc.in])
    AC_CONFIG_FILES([${pkg_config_install}/capdkrak.pc:src/capd/krak/capdkrak-install.pc.in])
   ]
)
# generate Makefile for krak, even if we do not use it. It is important for 'distclean' target in parent directory
AC_CONFIG_FILES([src/capd/krak/Makefile])


AC_SUBST([CAPD_MODULES])

AC_SUBST(capdAux_CFLAGS)
AC_SUBST(capdAux_LIBS)
AC_CONFIG_FILES([${pkg_config_dev}/capdAux.pc:src/capd/capdAux-dev.pc.in])
AC_CONFIG_FILES([${pkg_config_install}/capdAux.pc:src/capd/capdAux-install.pc.in])

AC_CONFIG_HEADERS([include/capd/config-capdAux.h])


AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([include/Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([src/capd/Makefile])
AC_CONFIG_FILES([src/capd/auxil/Makefile])

AC_CONFIG_SUBDIRS([examples])
AC_CONFIG_SUBDIRS([tests])


AC_OUTPUT
