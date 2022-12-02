include ${capdMake}/make/common_programs_makefile.mkf

noinst_PROGRAMS=alglibdemo
alglibdemo_SOURCES=alglibdemo.cpp
alglibdemo_libs=${capdvectalg_LIBS} ${capdmatrixAlgorithms_LIBS} ${capdAux_LIBS} ${capdExt_LIBS}
alglibdemo_LDADD=${alglibdemo_libs}
alglibdemo_DEPENDENCIES=${filter %.a %.la,${alglibdemo_libs}}
alglibdemo_CXXFLAGS= ${AM_CXXFLAGS} ${capdvectalg_CFLAGS} ${capdmatrixAlgorithms_CFLAGS} ${capdAux_CFLAGS} ${capdExt_CFLAGS}

