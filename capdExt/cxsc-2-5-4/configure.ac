AC_INIT([cxsc], [2.5.3])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_HEADERS([cxscconf.h])
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST
AC_PROG_INSTALL
AM_INIT_AUTOMAKE([-Wall -Wno-override foreign subdir-objects
	tar-pax no-dist-gzip dist-xz])
AC_PROG_CC
AM_PROG_CC_C_O
AM_PROG_AS
AC_PROG_CXX
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
LT_INIT([disable-static])

AH_TOP([
#ifndef _CXSCCONF_H_INCLUDED
#define _CXSCCONF_H_INCLUDED
])
AH_TEMPLATE([IBM_LINUX_C], [Compiling on a sane platform with a sane compiler])
AH_TEMPLATE([IBM_RS_6000_C], [Compiling with AIX compiler])
AH_TEMPLATE([SUN4_GNU_C], [Compiling with GCC under Solaris])
AH_TEMPLATE([SUN4_FORTE], [Compiling with Sun Forte (Solaris)])
AH_TEMPLATE([SUN_STUDIO_10], [Compiling with Sun Forte >= 5.7 (Solaris)])
AH_TEMPLATE([DEC_ALPHA_C], [Linux-ALPHA system])
AH_TEMPLATE([ROUND_C99_SAVE], [])
AH_TEMPLATE([ROUND_C99_QUICK], [])
AH_TEMPLATE([ROUND_C96_SAVE], [])
AH_TEMPLATE([ROUND_C96_QUICK], [])
AH_TEMPLATE([ROUND_ASM], [Employ ASM])
AH_TEMPLATE([CXSC_x86], [Use x86 ASM code if ROUND_ASM])
AH_TEMPLATE([CXSC_x86_64], [Use amd64 ASM code if ROUND_ASM])
AH_TEMPLATE([GNU_X86_64], [64-bit longs])
AH_BOTTOM([#endif])

AC_DEFINE([IBM_LINUX_C], [0])
AC_DEFINE([IBM_RS_6000_C], [0])
AC_DEFINE([SUN4_GNU_C], [0])
AC_DEFINE([SUN4_FORTE], [0])
AC_DEFINE([SUN_STUDIO_10], [0])
AC_DEFINE([DEC_ALPHA_C], [0])
AC_DEFINE([ROUND_C99_SAVE], [0])
AC_DEFINE([ROUND_C99_QUICK], [0])
AC_DEFINE([ROUND_C96_SAVE], [0])
AC_DEFINE([ROUND_C96_QUICK], [0])
AC_DEFINE([ROUND_ASM], [0])
AC_DEFINE([CXSC_x86], [0])
AC_DEFINE([CXSC_x86_64], [0])
AC_DEFINE([GNU_X86_64], [0])

AC_ARG_ENABLE([rounding],
	[
		AS_HELP_STRING([--enable-rounding=asm], [for hardware support for IEEE 754 arithmetic implemented in Assembler])
		AS_HELP_STRING([--enable-rounding=hardquick], [for hardware support for IEEE 754 arithmetic: for the slow and uncritical hardware IEEE arithmetic])
		AS_HELP_STRING([--enable-rounding=hardsafe], [for hardware support for IEEE 754 arithmetic: for the fast but critical hardware IEEE arithmetic])
		AS_HELP_STRING([--enable-rounding=soft], [for software emulations for directed rounded floating-point operations])
	],
	[rounding="$enableval"; echo "Rounding mode set to $enableval"],
	[rounding="hardsafe"]
)

# So make sure you pass in --host= if you want other bitness...
case "$host_cpu" in
x86_64|powerpc64|sparc64|ia64)
	# This macro's name is misleading; it has to do with sizeof(long)
	AC_DEFINE([GNU_X86_64], [1])
	BIT=64
	;;
*)
	BIT=32
	;;
esac

#
# Rounding choice
#
case "$host_os-$rounding" in
cygwin-asm)
	rounding="hardsafe";;
esac

case "$rounding-$host_cpu" in
asm-i[3456]86)
	AC_DEFINE([ROUND_ASM], [1])
	AC_DEFINE([CXSC_x86], [1])
	;;
asm-x86_64)
	AC_DEFINE([ROUND_ASM], [1])
	AC_DEFINE([CXSC_x86_64], [1])
	;;
hardquick-*)
	AC_DEFINE([ROUND_C99_QUICK], [1])
	;;
hardsafe-*)
	AC_DEFINE([ROUND_C99_SAVE], [1])
	;;
esac

#
# Transfer rounding choice to source file selection
#
AM_CONDITIONAL([HOST_DARWIN], [test "$host_os" = "darwin"])
AM_CONDITIONAL([USE_ASM_X86], [false])
AM_CONDITIONAL([USE_ASM_X86_64], [false])
case "$rounding-$host_cpu" in
asm-i[3456]86)
	AM_CONDITIONAL([USE_ASM_X86], [true])
	;;
asm-x86_64)
	AM_CONDITIONAL([USE_ASM_X86_64], [true])
	;;
esac

#
# Special SunOS deals
#
case "$host_os-$rounding-$CC" in
sunos-hardquick-gcc)
	AC_DEFINE([ROUND_C99_QUICK], [0])
	AC_DEFINE([ROUND_C96_QUICK], [1])
	;;
sunos-hardsoft-gcc)
	AC_DEFINE([ROUND_C99_SAVE], [0])
	AC_DEFINE([ROUND_C96_SAVE], [1])
	;;
esac
case "$host_os-$CC" in
sunos-cc)
	LIBS="$LIBS -lm9x";;
esac

#
# Compiler choice
#
case "$host_os-$host_cpu" in
linux-alpha)
	AC_DEFINE([DEC_ALPHA_C], [1])
	;;
linux-*|cygwin-*|darwin-i[3456]86|darwin-x86_64)
	AC_DEFINE([IBM_LINUX_C], [1])
	;;
sunos-*)
	AC_DEFINE([SUN4_GNU_C], [1])
	;;
darwin-*)
	AC_DEFINE([SUN4_GNU_C], [1])
	;;
aix-*)
	AC_DEFINE([IBM_RS_6000_C], [1])
	;;
esac
case "$host_os-$CC" in
sunos-cc)
	AC_DEFINE([SUN4_GNU_C], [0])
	AC_DEFINE([SUN4_FORTE], [1])
	dnl $SUN_FORTE_VERSION >= 5.7 then
	dnl AC_DEFINE([SUN4_FORTE], [0])
	dnl AC_DEFINE([SUN_STUDIO_10], [1])
	dnl uff!
	;;
esac

AM_CFLAGS="-Wall"
# our code is too shoddy after all
AS_IF([test "$CC" = "gcc"], [AM_CFLAGS="$AM_CFLAGS -fno-strict-aliasing"])
AM_CXXFLAGS="$AM_CFLAGS"
AC_SUBST([AM_CFLAGS])
AC_SUBST([AM_CXXFLAGS])
AC_CONFIG_FILES([Makefile CToolbox/Modules/Makefile
	CToolbox/Programs/Makefile CToolbox/test/Makefile benchmark/Makefile
	examples/Makefile src/Makefile src/rts/Makefile src/fi_lib/Makefile])
AC_OUTPUT

echo "Bitness  : $BIT"
echo "Rounding : $rounding"
echo "OS/CPU   : $host_os-$host_cpu"
