include ${capdMake}/make/common_programs_makefile.mkf

noinst_PROGRAMS=configFileReader configFileReader2

configFileReader_SOURCES=configFileReader.cpp
configFileReader_libs=${capdauxil_LIBS}
configFileReader_LDADD=${configFileReader_libs}
configFileReader_DEPENDENCIES=${filter %.a %.la,${configFileReader_libs}}
configFileReader_CXXFLAGS= ${AM_CXXFLAGS} ${capdauxil_CFLAGS}

configFileReader2_SOURCES=configFileReader2.cpp
configFileReader2_libs=${capdauxil_LIBS}
configFileReader2_LDADD=${configFileReader2_libs}
configFileReader2_DEPENDENCIES=${filter %.a %.la,${configFileReader2_libs}}
configFileReader2_CXXFLAGS= ${AM_CXXFLAGS} ${capdauxil_CFLAGS}

