include ${capdMake}/make/common_programs_makefile.mkf

check_PROGRAMS=divtst intvtst
#ioTest # this test is broken for a long time
TESTS=${check_PROGRAMS}

#XFAIL_TESTS=ioTest

divtst_SOURCES=divtst.cpp
divtst_libs=${capdAlg_LIBS} ${capdExt_LIBS}
divtst_LDADD=${divtst_libs}
divtst_DEPENDENCIES=${filter %.a %.la,${divtst_libs}}
divtst_CXXFLAGS= ${AM_CXXFLAGS} ${capdintervals_CFLAGS} ${capdExt_CFLAGS} ${capdExt_CFLAGS}

intvtst_SOURCES=intvtst.cpp
intvtst_libs=${capdAlg_LIBS} ${capdExt_LIBS}
intvtst_LDADD=${intvtst_libs}
intvtst_DEPENDENCIES=${filter %.a %.la,${intvtst_libs}}
intvtst_CXXFLAGS= ${AM_CXXFLAGS} ${capdintervals_CFLAGS} ${capdExt_CFLAGS} ${capdExt_CFLAGS}

# ioTest_SOURCES=ioTest.cpp
# ioTest_libs=${capdAlg_LIBS} ${capdExt_LIBS}
# ioTest_LDADD=${ioTest_libs}
# ioTest_DEPENDENCIES=${filter %.a %.la,${ioTest_libs}}
# ioTest_CXXFLAGS=${capdintervals_CFLAGS} ${capdExt_CFLAGS}
