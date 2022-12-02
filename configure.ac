#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.


#AC_PREREQ([2.61])
m4_define([capdMake_PATH],[capdMake])
#m4_esyscmd(capdMake_PATH[/utils/capd_library_version.sh] capdMake_PATH[/../] keep))

m4_include(capdMake_PATH[/make/capd_version.m4])

AC_INIT([capd],[CAPD_VERSION_NUMBER],[http://capd.ii.uj.edu.pl],[],[http://capd.ii.uj.edu.pl])

AC_CONFIG_MACRO_DIRS([capdMake/ac_aux/m4 .]) # dot . for old autotools in Debian
AC_CONFIG_AUX_DIR(capdMake/ac_aux)
default_top_capd=${srcdir}
top_capd_build=.

m4_include([capdMake/make/configure_common.ac])

AM_INIT_AUTOMAKE(capd_automake_init_options)


if test "$cache_file" = /dev/null; then
  cache_file="${top_capd_build}/capd.cache"
  if test -e $cache_file; then
     ${RM} $cache_file
   fi
fi


AC_CONFIG_FILES([Makefile])

CAPD_PACKAGES=""

if test -e "${srcdir}/capdExt"; then
  AC_CONFIG_SUBDIRS([capdExt])
  CAPD_PACKAGES="$CAPD_PACKAGES  capdExt "
fi

if test -e "${srcdir}/capdAux"; then
  AC_CONFIG_SUBDIRS([capdAux])
  CAPD_PACKAGES="$CAPD_PACKAGES  capdAux "
fi

if test -e "${srcdir}/capdAlg"; then
  AC_CONFIG_SUBDIRS([capdAlg])
  CAPD_PACKAGES="$CAPD_PACKAGES  capdAlg "
fi

if test -e "${srcdir}/capdExtHom"; then
  AC_CONFIG_SUBDIRS([capdExtHom])
  CAPD_PACKAGES="$CAPD_PACKAGES  capdExtHom "
fi

if test -e "${srcdir}/capdRedHom"; then
  AC_CONFIG_SUBDIRS([capdRedHom])
  CAPD_PACKAGES="$CAPD_PACKAGES  capdRedHom "
fi

if test -e "${srcdir}/capdDynSys"; then
  AC_CONFIG_SUBDIRS([capdDynSys])
  CAPD_PACKAGES="$CAPD_PACKAGES  capdDynSys"
else
  if test -e "${srcdir}/capdDynSys4"; then
    AC_CONFIG_SUBDIRS([capdDynSys4])
    CAPD_PACKAGES="$CAPD_PACKAGES  capdDynSys4"
  fi
fi

AC_CONFIG_SUBDIRS([capdMake/libcapd])

AC_SUBST([CAPD_PACKAGES])


AC_PATH_PROG([DEBUILD], [dpkg-buildpackage], [notfound])
AC_PATH_PROG([TEST], [test])
AC_PATH_PROG([TAR], [tar gtar])   # tar
AC_PATH_PROG([XARGS], [xargs])    # xargs
AC_PATH_PROG([CHMOD], [chmod])    # chmod
AC_PATH_PROG([CAT], [cat])        # cat
AC_PATH_PROG([ECHO], [echo])      # echo
AC_PATH_PROG([CP], [cp])          # cp
AC_PATH_PROG([DATE], [date])      # date
AC_PATH_PROG([MV], [mv])          # mv
AC_PATH_PROG([TAIL], [tail])          # tail


AC_OUTPUT
