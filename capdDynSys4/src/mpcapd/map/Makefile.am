include ${capdMake}/make/common_makefile.mkf

noinst_LTLIBRARIES=libmpcapdmap.la

libmpcapdmap_la_SOURCES=  Function.cpp  Map.cpp 
libmpcapdmap_la_CXXFLAGS= ${AM_CXXFLAGS}  ${mpcapdmap_CFLAGS} ${mpcapdAlg_CFLAGS} 

pkgconfigdir = @PKG_CONFIG_DIR@
pkgconfig_DATA =
