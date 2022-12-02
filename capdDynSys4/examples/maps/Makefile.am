include ${capdMake}/make/common_programs_makefile.mkf

noinst_PROGRAMS=mapExample hessianExample jetExample jetIndexingExample jetTransportExample

mapExample_SOURCES=mapExample.cpp
mapExample_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS}
mapExample_libs=${capdDynSys_LIBS}
mapExample_LDADD=${mapExample_libs}
mapExample_DEPENDENCIES=${filter %.a %.la,${mapExample_libs}}

hessianExample_SOURCES=hessianExample.cpp
hessianExample_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS}
hessianExample_libs=${capdDynSys_LIBS}
hessianExample_LDADD=${hessianExample_libs}
hessianExample_DEPENDENCIES=${filter %.a %.la,${hessianExample_libs}}

jetExample_SOURCES=jetExample.cpp
jetExample_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS}
jetExample_libs=${capdDynSys_LIBS}
jetExample_LDADD=${jetExample_libs}
jetExample_DEPENDENCIES=${filter %.a %.la,${jetExample_libs}}

jetIndexingExample_SOURCES=jetIndexingExample.cpp
jetIndexingExample_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS}
jetIndexingExample_libs=${capdDynSys_LIBS}
jetIndexingExample_LDADD=${jetIndexingExample_libs}
jetIndexingExample_DEPENDENCIES=${filter %.a %.la,${jetIndexingExample_libs}}

jetTransportExample_SOURCES=jetTransportExample.cpp
jetTransportExample_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS}
jetTransportExample_libs=${capdDynSys_LIBS}
jetTransportExample_LDADD=${jetExample_libs}
jetTransportExample_DEPENDENCIES=${filter %.a %.la,${jetTransportExample_libs}}

EXTRA_DIST=maps.dox
