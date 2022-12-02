include ${capdMake}/make/common_programs_makefile.mkf

noinst_PROGRAMS=DPoincareMapExample DPoincareMapDerivativeExample IPoincareMapExample IPoincareMapDerivativeExample

DPoincareMapExample_SOURCES=DPoincareMapExample.cpp
DPoincareMapExample_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS}
DPoincareMapExample_libs=${capdDynSys_LIBS}
DPoincareMapExample_LDADD=${DPoincareMapExample_libs}
DPoincareMapExample_DEPENDENCIES=${filter %.a %.la,${DPoincareMapExample_libs}}

DPoincareMapDerivativeExample_SOURCES=DPoincareMapDerivativeExample.cpp
DPoincareMapDerivativeExample_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS}
DPoincareMapDerivativeExample_libs=${capdDynSys_LIBS}
DPoincareMapDerivativeExample_LDADD=${DPoincareMapDerivativeExample_libs}
DPoincareMapDerivativeExample_DEPENDENCIES=${DPoincareMapDerivativeExample_libs:-%=}

IPoincareMapExample_SOURCES=IPoincareMapExample.cpp
IPoincareMapExample_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS}
IPoincareMapExample_libs=${capdDynSys_LIBS}
IPoincareMapExample_LDADD=${IPoincareMapExample_libs}
IPoincareMapExample_DEPENDENCIES=${IPoincareMapExample_libs:-%=}

IPoincareMapDerivativeExample_SOURCES=IPoincareMapDerivativeExample.cpp
IPoincareMapDerivativeExample_CXXFLAGS= ${AM_CXXFLAGS} ${capdDynSys_CFLAGS}
IPoincareMapDerivativeExample_libs=${capdDynSys_LIBS}
IPoincareMapDerivativeExample_LDADD=${IPoincareMapDerivativeExample_libs}
IPoincareMapDerivativeExample_DEPENDENCIES=${IPoincareMapDerivativeExample_libs:-%=}

EXTRA_DIST=poincare-nonrigorous.dox poincare-rigorous.dox poincare.dox 
