AC_INIT([filib],[3.0.2])
AC_CONFIG_HEADERS(config.h)
AM_INIT_AUTOMAKE
AC_PROG_CC
CXXFLAGS_SAVE=" -I${PWD}/${srcdir} ${CXXFLAGS}"
CXXFLAGS=
AC_PROG_CXX
CXXFLAGS=${CXXFLAGS_SAVE}
AC_PROG_RANLIB
#AM_DISABLE_SHARED
#AC_PROG_LIBTOOL
LT_INIT
AC_HEADER_STDC
#m4_warn([obsolete],
#[The preprocessor macro `STDC_HEADERS' is obsolete.
#  Except in unusual embedded environments, you can safely include all
#  ISO C90 headers unconditionally.])dnl
# Autoupdate added the next two lines to ensure that your configure
# script's behavior did not change.  They are probably safe to remove.
#AC_CHECK_INCLUDES_DEFAULT
AC_PROG_EGREP

AC_CHECK_HEADERS([sys/types.h])

AC_DEFINE(FILIB_EXTENDED, [], [FILIB_EXTENDED])
# AC_DEFINE(FILIB_PRED_SUCC_TABLES, [], [Use talbes to speedup the computation of the pred and succ function.])

fparchitecture=

AC_ARG_WITH(
	[force-arch],
	[  --with-force-arch=[[no]]
                          override automatic floating point architecture selection [[default=no]]
                          available options are x87,SSE],
        [if [[ ! -z "${withval}" ]] ; then fparchitecture=${withval} ; fi ]
)

if [[ -z "${fparchitecture}" ]] ; then
AC_MSG_CHECKING([for floating point architecture])
	case "$host" in
	  x86_64-*-*)
	    fparchitecture=SSE
	    ;;
	  i?86-*-*)
	    fparchitecture=x87
	    ;;
	  *)
	    fparchitecture="unknown floating point architecture in host \"$host\""
	    ;;
	esac
AC_MSG_RESULT([$fparchitecture])
else
	echo "User forced floating point architecture to $fparchitecture."

	case "$fparchitecture" in
		SSE)
			case "$host" in
			  x86_64-*-*)
			    ;;
			  i686-*-*)
			    CXXFLAGS="${CXXFLAGS} -msse -march=pentium4"
			    ;;
			  *)
			    echo "WARNING: forced floating point architecture to SSE on unknown platform."
			    ;;
			esac
			;;
		x87)
			;;
		*)
			;;
	esac
fi

define_have_sse="// #define HAVE_SSE"
define_have_x87="// #define HAVE_X87"

case "${fparchitecture}" in
  SSE)
    define_have_sse="#define HAVE_SSE"
    echo "Using SSE floating point architecture."
    WANT_ROUNDING_MATH=yes
    WANT_FLOAT_STORE=no
    ;;
  x87)
    define_have_x87="#define HAVE_X87"
    echo "Using x87 floating point architecture."
    WANT_ROUNDING_MATH=yes
    WANT_FLOAT_STORE=yes
    ;;
  *)
    echo "Unknown floating point architecture."
    WANT_ROUNDING_MATH=yes
    WANT_FLOAT_STORE=yes
    ;;
esac

if [[ "$GXX" = "yes" ]] ; then
	AC_LANG_PUSH([C++])

        CXXFLAGS_SAVE="${CXXFLAGS}"
	CXXFLAGS="${CXXFLAGS} -frounding-math"
	AC_MSG_CHECKING([whether gcc accepts -frounding-math])
	AC_LINK_IFELSE([AC_LANG_SOURCE([int main() {}])],
	have_frounding_math=yes,
	have_frounding_math=no
	)
	AC_MSG_RESULT($have_frounding_math)
	if [[ "$have_frounding_math" = "no" -o "$WANT_ROUNDING_MATH" = "no" ]] ; then CXXFLAGS="$CXXFLAGS_SAVE" ; fi

        CXXFLAGS_SAVE="${CXXFLAGS}"
	CXXFLAGS="${CXXFLAGS} -ffloat-store"
	AC_MSG_CHECKING([whether gcc accepts -ffloat-store])
	AC_LINK_IFELSE([AC_LANG_SOURCE([int main() {}])],
	have_float_store=yes,
	have_float_store=no
	)
	AC_MSG_RESULT($have_float_store)
	if [[ "$have_float_store" = "no" -o "$WANT_FLOAT_STORE" = "no" ]] ; then CXXFLAGS="$CXXFLAGS_SAVE" ; fi

	AC_LANG_POP
fi

AC_SUBST([define_have_sse])
AC_SUBST([define_have_x87])

AC_LANG_PUSH([C++])
AC_MSG_CHECKING([whether the boost random library is installed])
  AC_LINK_IFELSE([AC_LANG_SOURCE([#include <boost/random.hpp>
int main() {

}])],
    have_boost=yes,
    have_boost=no
  )
  AC_MSG_RESULT($have_boost)
AC_LANG_POP

benchmark_subdir=""

if [[ "$have_boost" = "yes" ]] ; then
	benchmark_subdir="benchmark"
	CPPFLAGS="${CPPFLAGS} -DHAVE_BOOST"
fi

AC_SUBST([benchmark_subdir])

AC_CONFIG_FILES([Makefile libprim/Makefile rounding_control/rounding_control_config.hpp])
AC_OUTPUT
# examples/Makefile test/Makefile benchmark/Makefile

echo
echo "Note on using the library:"
echo
echo "Please set your CXXFLAGS to \"${CXXFLAGS}\" when compiling programs using the installed filib++."
echo
