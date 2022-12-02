include ${capdMake}/make/common_programs_makefile.mkf

noinst_PROGRAMS=mpExample1 mpExample2 mpExample3 mpExample4 multiprecExample


mpExample1_SOURCES=mpExample1.cpp
mpExample1_libs=${mpcapddynsys_LIBS} ${mpfr_LIBS}
mpExample1_LDADD=${mpExample1_libs}
mpExample1_DEPENDENCIES=${filter %.a %.la,${mpExample1_libs}}
mpExample1_CXXFLAGS= ${AM_CXXFLAGS} ${mpcapdDynSys_CFLAGS} ${mpfr_CFLAGS}

mpExample2_SOURCES=mpExample2.cpp
mpExample2_CXXFLAGS= ${AM_CXXFLAGS} ${mpcapdDynSys_CFLAGS} ${mpfr_CFLAGS}
mpExample2_libs=${mpcapdDynSys_LIBS} ${mpfr_LIBS}
mpExample2_LDADD=${mpExample2_libs}
mpExample2_DEPENDENCIES=${filter %.a %.la,${mpExample2_libs}}

mpExample3_SOURCES=mpExample3.cpp
mpExample3_CXXFLAGS= ${AM_CXXFLAGS} ${mpcapdDynSys_CFLAGS} ${mpfr_CFLAGS}
mpExample3_libs=${mpcapdDynSys_LIBS}  ${mpfr_LIBS}
mpExample3_LDADD=${mpExample3_libs}
mpExample3_DEPENDENCIES=${filter %.a %.la,${mpExample3_libs}}

mpExample4_SOURCES=mpExample4.cpp
mpExample4_CXXFLAGS= ${AM_CXXFLAGS} ${mpcapdDynSys_CFLAGS}  ${mpfr_CFLAGS}
mpExample4_libs=${mpcapdDynSys_LIBS} ${mpfr_LIBS}
mpExample4_LDADD=${mpExample4_libs}
mpExample4_DEPENDENCIES=${filter %.a %.la,${mpExample4_libs}}

# mpExample6_SOURCES=mpExample6.cpp
# mpExample6_CXXFLAGS= ${AM_CXXFLAGS} ${mpcapdDynSys_CFLAGS}  ${mpfr_CFLAGS}
# mpExample6_libs=${mpcapdDynSys_LIBS} ${mpfr_LIBS}
# mpExample6_LDADD=${mpExample6_libs}
# mpExample6_DEPENDENCIES=${filter %.a %.la,${mpExample6_libs}}



# mpExample5_SOURCES=mpExample5.cpp
# mpExample5_CXXFLAGS= ${AM_CXXFLAGS} ${mpcapdDynSys_CFLAGS}  ${mpfr_CFLAGS}
# mpExample5_libs=${mpcapdDynSys_LIBS} ${mpfr_LIBS}
# mpExample5_LDADD=${mpExample5_libs}
# mpExample5_DEPENDENCIES=${filter %.a %.la,${mpExample5_libs}}

multiprecExample_SOURCES=multiprecExample.cpp
multiprecExample_CXXFLAGS= ${AM_CXXFLAGS} ${mpcapdDynSys_CFLAGS} ${mpfr_CFLAGS}
multiprecExample_libs=${mpcapdDynSys_LIBS}  ${mpfr_LIBS}
multiprecExample_LDADD=${multiprecExample_libs}
multiprecExample_DEPENDENCIES=${filter %.a %.la,${multiprecExample_libs}}

EXTRA_DIST=multiprec.dox
