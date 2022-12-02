include ${capdMake}/make/common_makefile.mkf
ACLOCAL_AMFLAGS= -I ../capdMake/ac_aux/m4

SUBDIRS=

if WITH_FILIB_INTERNAL
	SUBDIRS+=filibsrc
endif

SUBDIRS+=src include


.PHONY: lib

lib:
	@for dir in ${SUBDIRS}; do \
	   echo "Making lib in $$dir"; \
	   $(MAKE) $(AM_MAKEFLAGS) -C $$dir ; \
	done;
