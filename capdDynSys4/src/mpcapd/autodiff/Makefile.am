include ${capdMake}/make/common_makefile.mkf

noinst_LTLIBRARIES=libmpcapdautodiff.la
libmpcapdautodiff_la_SOURCES= DagIndexer.cpp
libmpcapdautodiff_la_LDFLAGS=${LIBCAPD_LDFLAGS}
libmpcapdautodiff_la_CXXFLAGS= ${AM_CXXFLAGS}  ${mpcapdautodiff_CFLAGS} ${mpcapdAlg_CFLAGS} ${capdAux_CFLAGS}
