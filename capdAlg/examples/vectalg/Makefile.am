include ${capdMake}/make/common_programs_makefile.mkf

noinst_PROGRAMS=matrixExample  normExample  vectorExample

matrixExample_SOURCES=matrixExample.cpp
matrixExample_libs=${capdintervals_LIBS} ${capdExt_LIBS}
matrixExample_LDADD=${matrixExample_libs}
matrixExample_DEPENDENCIES=${filter %.a %.la,${matrixExample_libs}}
matrixExample_CXXFLAGS= ${AM_CXXFLAGS} ${capdintervals_CFLAGS} ${capdExt_CFLAGS}

normExample_SOURCES=normExample.cpp
normExample_libs=${capdvectalg_LIBS} ${capdintervals_LIBS} ${capdExt_LIBS}
normExample_LDADD=${normExample_libs}
normExample_DEPENDENCIES=${filter %.a %.la,${normExample_libs}}
normExample_CXXFLAGS= ${AM_CXXFLAGS} ${capdvectalg_CFLAGS} ${capdintervals_CFLAGS} ${capdExt_CFLAGS}

vectorExample_SOURCES=vectorExample.cpp
vectorExample_libs=${capdvectalg_LIBS} ${capdintervals_LIBS} ${capdExt_LIBS}
vectorExample_LDADD=${vectorExample_libs}
vectorExample_DEPENDENCIES=${filter %.a %.la,${vectorExample_libs}}
vectorExample_CXXFLAGS= ${AM_CXXFLAGS} ${capdvectalg_CFLAGS} ${capdintervals_CFLAGS} ${capdExt_CFLAGS}


