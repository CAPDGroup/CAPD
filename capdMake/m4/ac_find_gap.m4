#   LINBOXING - ac_find_gap.m4
#   Paul Smith
#
#   Copyright (C)  2007-2008
#   National University of Ireland Galway
#   Copyright (C)  2011
#   University of St Andrews
#
#   This file is part of the linboxing GAP package.
#
#   The linboxing package is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or (at your
#   option) any later version.
#
#   The linboxing package is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
#   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
#   more details.
#
#   You should have received a copy of the GNU General Public License along with
#   this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

# Find the location of GAP
# Sets GAPROOT, GAPARCH and GAP_CPPFLAGS appropriately
#####################################################

AC_DEFUN([AC_FIND_GAP],
[
  AC_LANG_PUSH(C)

  # Make sure CDPATH is portably set to a sensible value
  CDPATH=${ZSH_VERSION+.}:

  GAP_CPPFLAGS=""

  ######################################
  # Find the GAP root directory by
  # checking for the sysinfo.gap file
  AC_MSG_CHECKING([for GAP root directory])

  #Allow the user to specify the location of GAP
  #
  AC_ARG_WITH(gaproot,
    [AC_HELP_STRING([--with-gaproot=<path>], [specify root of GAP installation])],
    [GAPROOT=$withval])
  DEFAULT_GAPROOTS="../.. ${HOME}/work/gap ${HOME}/work/gap4r8 ${HOME}/work/gap4r7 ${HOME}/work/gap4r6 ${HOME}/gap ${HOME}/gap4r8 ${HOME}/gap4r7 ${HOME}/gap4r6 ${HOME}/gap4r5 ${HOME}/gap4r4 /usr/local/lib/gap /usr/lib/gap /usr/local/lib/gap4r8 /usr/lib/gap4r8 /usr/local/lib/gap4r7 /usr/lib/gap4r7 /usr/local/lib/gap4r6 /usr/lib/gap4r6 /usr/local/lib/gap4r5 /usr/lib/gap4r5 /usr/local/lib/gap4r4 /usr/lib/gap4r4"

  SYSINFO="sysinfo.gap"
  havesysinfo=0
  if test -n "$GAPROOT"; then
    # Make sure the path given is stored as an absolute path
    GAPROOT=`cd $GAPROOT > /dev/null 2>&1 && pwd`

    if test -e ${GAPROOT}/${SYSINFO}; then
      havesysinfo=1
    fi
  else
    AS_ECHO(["Looking for GAP in ${DEFAULT_GAPROOTS}"])
    # Otherwise try likely directories
    for GAPROOT in ${DEFAULT_GAPROOTS}
    do
      # Convert the path to absolute
      GAPROOT=`cd $GAPROOT > /dev/null 2>&1 && pwd`
      if test -e ${GAPROOT}/${SYSINFO}; then
        havesysinfo=1
        break
      fi
    done
  fi

  if test "x$havesysinfo" != "x1"; then
    echo ""
    echo "**********************************************************************"
    echo ""
    echo "  Cannot find your GAP installation. Please specify the location of"
    echo "  GAP's root directory using --with-gaproot=<path>"
    echo ""
    echo "**********************************************************************"
    echo ""

    AC_MSG_NOTICE([Unable to find GAP root directory])
    HAVE_GAP=no
  else
    HAVE_GAP=yes
    AC_MSG_RESULT([${GAPROOT}])

  #####################################
  # First check the version

  AC_REQUIRE([AC_PROG_GREP])
  AC_REQUIRE([AC_PROG_SED])

  GAPVERSION="Unknown"

  AC_MSG_CHECKING([for GAP version])
  GAPSYSTEM_G="$GAPROOT/lib/system.g"
  if test -r "$GAPSYSTEM_G"; then
    GAPVERSION=`${GREP} ' Version *:= *"' ${GAPSYSTEM_G} | ${SED} 's|Version *:= *"\(.*\)",.*|\1|'`
    AC_MSG_RESULT([${GAPVERSION}])

    AC_MSG_CHECKING([for GAP version >= 4.5])
    GAP45=no
    AX_COMPARE_VERSION([$GAPVERSION], [ge], [4.5], [
      GAP45=yes
    ], [
      AX_COMPARE_VERSION([$GAPVERSION], [eq], [4.dev], [GAP45=yes])])

    AC_MSG_RESULT($GAP45)

    if test "x$GAP45" = "xyes"; then
      AC_DEFINE(GAP_4_5, , "We have GAP version 4.5 (or higher)")
    fi
  else
    AC_MSG_RESULT([Not found])

    echo ""
    echo "**********************************************************************"
    echo "  ERROR"
    echo ""
    echo "  Found a GAP installation at $GAPROOT"
    echo "  but could find information about the GAP version in the file"
    echo "  ${GAPSYSTEM_G}"
    echo "  This file should be present: please check your GAP installation."
    echo ""
    echo "**********************************************************************"
    echo ""

    AC_MSG_WARN([Unable to find GAP system.g])
    HAVE_GAP=no
  fi

  if test "x$GAPARCH" = "xUnknown"; then
    echo ""
    echo "**********************************************************************"
    echo "  ERROR"
    echo ""
    echo "  Found a GAP installation at $GAPROOT"
    echo "  but could find information about GAP's architecture in the file"
    echo "  ${GAPROOT}/${SYSINFO}"
    echo "  This file should be present: please check your GAP installation."
    echo "**********************************************************************"
    echo ""

    AC_MSG_WARN([Unable to find GAParch information.])
    HAVE_GAP=no
  fi


  #####################################
  # Now find the architecture


  GAPARCH="Unknown"

  AC_MSG_CHECKING([for GAP architecture])
  GAPARCH=`${GREP} GAParch= ${GAPROOT}/${SYSINFO} | ${SED} 's|^GAParch=\(.*\)|\1|'`
  AC_MSG_RESULT([${GAPARCH}])

  if test "x$GAPARCH" = "xUnknown"; then
    echo ""
    echo "**********************************************************************"
    echo "  ERROR"
    echo ""
    echo "  Found a GAP installation at $GAPROOT"
    echo "  but could find information about GAP's architecture in the file"
    echo "  ${GAPROOT}/${SYSINFO}"
    echo "  This file should be present: please check your GAP installation."
    echo ""
    echo "**********************************************************************"
    echo ""

    AC_MSG_WARN([Unable to find GAParch information.])
    HAVE_GAP=no
  fi



  ######################################
  # Check for the GAP config.h

  OLD_CPPFLAGS="$CFLAGS"

  GAPCONFIGPATH="$GAPROOT/bin/$GAPARCH"
  CPPFLAGS="$CPPFLAGS -I$GAPCONFIGPATH"

  AC_CHECK_HEADERS([config.h],
    [
      GAP_CPPFLAGS="$GAP_CPPFLAGS -DCONFIG_H -I$GAPCONFIGPATH"
      HAVECONFIGH=true
    ],
    [
      HAVECONFIGH=false
    ])

  CPPFLAGS="$OLD_CPPFLAGS"

  #####################################
  # Now check for the GAP header files

  #Remember the input CPPFLAGS
  OLD_CPPFLAGS="$CPPFLAGS"

  AC_MSG_CHECKING([for GAP include files])
  CPPFLAGS="$OLD_CPPFLAGS -I$GAPROOT  -DCONFIG_H -I$GAPCONFIGPATH"

  AC_COMPILE_IFELSE(
    [AC_LANG_SOURCE([[#include "src/compiled.h"
      Obj Func(Obj Self){return (Obj)0;}]])],
    [a=1])

  if test "x$a" = "x1"; then
    GAP_CPPFLAGS="-I$GAPROOT"
    AC_MSG_RESULT([$GAPROOT/src])
  else
    AC_MSG_RESULT([Not found])

    echo ""
    echo "**********************************************************************"
    echo "  ERROR"
    echo ""
    echo "  Failed to find the GAP source header files in the src/ subdirectory"
    echo "  of your GAP root directory. The expected location was"
    echo "  $GAPROOT/src"
    echo ""
    echo "  The linboxing kernel build process expects to find the normal GAP "
    echo "  root directory structure as it is after building GAP itself, and in "
    echo "  particular the files <gaproot>/sysinfo.gap, <gaproot>/src/<includes>"
    echo "  and <gaproot>/bin/<architecture>/bin/config.h. Please make sure that"
    echo "  your GAP root directory structure conforms to this layout, or give"
    echo "  the correct GAP root using  --with-gaproot=<path>"
    echo "**********************************************************************"
    echo ""

    AC_MSG_WARN([Unable to find GAP include files in /src subdirectory])
    HAVE_GAP=no
  fi

  #Reset CPPFLAGS
  CPPFLAGS="$OLD_CPPFLAGS"


  ######################################
  # Try to find out whether GAP was compiled as 32-bit or 64-bit

  AC_MSG_CHECKING([sysinfo.gap for 32-bit or 64-bit GAP])
  # If we are GAP 4.5 or above it will be stored in sysinfo.gap
  BUILDMODE=`${GREP} GAParch_abi= ${GAPROOT}/${SYSINFO} | ${SED} 's|^GAParch_abi=\(..\)-bit|\1|'`

  if test "x$BUILDMODE" = "x"; then
    AC_MSG_RESULT([not found])
    if test "x$HAVECONFIGH" = "xtrue"; then
      # otherwise ask config.h (if present) for the size of a void*
      AC_LANG_PUSH(C)
      OLD_CPPFLAGS=$CPPFLAGS
      CPPFLAGS="$OLDCPPFLAGS -I$GAPCONFIGPATH"
      AC_MSG_CHECKING([config.h for 32-bit or 64-bit GAP])
      AC_RUN_IFELSE(
      [
        AC_LANG_PROGRAM(
        [[#include <config.h> ]],
        [[
          if(SIZEOF_VOID_P == 8)
            return 0;
          else
            return -1;
        ]])
      ], [BUILDMODE=64],
      [
        AC_RUN_IFELSE(
        [
          AC_LANG_PROGRAM(
          [[#include <config.h> ]],
          [[
            if(SIZEOF_VOID_P == 4)
              return 0;
            else
              return -1;
          ]])
        ], [BUILDMODE=32])
      ])
      AC_LANG_POP(C)
      CPPFLAGS=$OLD_CPPFLAGS
    fi

    if test "x$BUILDMODE" = "x"; then
      AC_MSG_RESULT([not found])
      AC_CHECK_SIZEOF([void *],4)
      bits64=""
      AC_MSG_CHECKING([for 64-bit machine])
      if test "x${ac_cv_sizeof_void_p}" = "x8"; then
        BUILDMODE=64
        AC_MSG_RESULT([yes])
        bits64="a 64-bit machine"
      else
        AC_MSG_RESULT([no])
        bits64="not a 64-bit machine"
      fi
      echo ""
      echo "----------------------------------------------------------------------"
      echo "  WARNING"
      echo ""
      echo "  Failed to find work out whether GAP thinks the system is 32- or"
      echo "  64-bit. As a result, the linboxing configure script has done its"
      echo "  own test and thinks it is $bits64. You can proceed to build the"
      echo "  kernel module, but it is important that you use"
      echo ""
      echo "  gap> TestLinboxing()"
      echo ""
      echo "  to make sure that linboxing and GAP agree about the size of small "
      echo "  integers."
      echo ""
      echo "  To avoid this warning, please check your GAP installation, and"
      echo "  that the architecture reported in <gaproot>/sysinfo.gap agrees"
      echo "  with the subdirectory of <gaproot>/bin"
      echo "----------------------------------------------------------------------"
      echo ""
    else
      AC_MSG_RESULT([$BUILDMODE-bit])
    fi
  else
    AC_MSG_RESULT([$BUILDMODE-bit])
  fi

  if test "x$BUILDMODE" = "x64"; then
    GAP_CPPFLAGS="$GAP_CPPFLAGS -DSYS_IS_64_BIT -m64"
    AC_DEFINE(SIZEOF_VOID_P, 8, [The size of `void *' as reported by GAP])
  elif test "x$BUILDMODE" = "x32"; then
    AC_DEFINE(SIZEOF_VOID_P, 4, [The size of `void *' as reported by GAP])
    GAP_CPPFLAGS="$GAP_CPPFLAGS -m32"
  else
    AC_MSG_WARN([Unrecognised BUILDMODE: $BUILDMODE. This is probably an error
      in the configure script. Please contact the package maintainer.])
    HAVE_GAP=no
  fi

  ######################################
  # Try to find out whether GAP was compiled with GMP

  AC_MSG_CHECKING([whether GAP was built to use GMP])

  # We have to check the Makefile
  WITHGMP=`${GREP} USE_GMP ${GAPROOT}/Makefile`

  if test "x$WITHGMP" = "x"; then
    AC_MSG_RESULT([no])
  else
    AC_MSG_RESULT([yes])
    AC_DEFINE(USE_GMP, , [GAP was built with GMP])
  fi

  AC_SUBST(BUILDMODE)
  AC_SUBST(GAPARCH)
  AC_SUBST(GAPROOT)
  AC_SUBST(GAP_CPPFLAGS)


  AC_LANG_POP(C)

  fi #   if test "x$havesysinfo" != "x1"; then

  if test "$HAVE_GAP" == "yes"; then
     HAVE_GAP=true
     GAP_CFLAGS="-I${GAPROOT}"
     AC_SUBST([GAP_CFLAGS])
  else
     HAVE_GAP=false
  fi

  AM_CONDITIONAL([HAVE_GAP], [$HAVE_GAP])

])
