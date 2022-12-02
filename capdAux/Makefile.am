include ${capdMake}/make/common_makefile.mkf
ACLOCAL_AMFLAGS= -I ../capdMake/ac_aux/m4

SUBDIRS=include src tests

if WITH_CAPD_EXAMPLES
  SUBDIRS+= examples
endif

.PHONY: lib

lib:
	$(MAKE) $(AM_MAKEFLAGS) -C src
