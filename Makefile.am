include ${capdMake}/make/common_makefile.mkf
ACLOCAL_AMFLAGS= -I capdMake/ac_aux/m4

SUBDIRS=. ${CAPD_PACKAGES} capdMake/libcapd

.PHONY: lib deb release

distclean-local:
	-rm -fr lib capd.cache

lib:
	@for dir in ${CAPD_PACKAGES}; do \
	   echo "Making lib in $$dir"; \
	   $(MAKE) $(AM_MAKEFLAGS) -C $$dir lib ; \
	   done;
	@$(MAKE) $(AM_MAKEFLAGS) -C "capdMake/libcapd" lib;



DEBSOURCEPKG=$(PACKAGE)_$(VERSION).orig.tar.gz
DEBSOURCEDIR=$(PACKAGE)-$(VERSION)
DEBFILES=$(top_srcdir)/capdMake/debian

deb: dist
	@$(ECHO) "Building Debian package..."; \
	if $(TEST) ! -d $(DEBFILES); then \
		$(ECHO) "The files needed for building a Debian package are not" \
		"included by default in the distribution. To build a package, check" \
		"out the project from source control."; \
		exit 1; \
	fi; \
	$(MV) $(PACKAGE)-$(VERSION).tar.gz $(DEBSOURCEPKG) && \
	$(TAR) --extract --gunzip --file $(DEBSOURCEPKG) && \
	$(MKDIR_P) $(DEBSOURCEDIR)/debian && \
	$(CP) -fr \
		$(DEBFILES) \
		$(DEBSOURCEDIR)/ && \
	cd $(DEBSOURCEDIR) && \
	$(DEBUILD) -j$(shell grep -c ^processor /proc/cpuinfo)  -rfakeroot -D -us -uc; \
	cd $(srcdir); \
	rm --force $(DEBSOURCEPKG); \
	rm --recursive --force $(DEBSOURCEDIR)



EXTRA_DIST = capdMake/make/capd_version.m4 capdMake/make/common_makefile.mkf capdMake/make/common_programs_makefile.mkf capdMake/make/configure_common.ac capdMake/make/configure_common_boost.ac capdMake/make/configure_common_cxx11.ac capdMake/make/configure_common_mathematica.ac capdMake/make/configure_common_mpfr.ac capdMake/make/configure_common_mpi.ac capdMake/make/configure_common_pari.ac capdMake/make/configure_common_tbb.ac capdMake/make/configure_common_x.ac capdMake/debian capdMake/capd_library_version_number.m4 capdMake/capd_version_number.m4 capdMake/utils/capd_library_version.sh  README.md capdMake/docs capdAlg/docs capdMake/examples


dist-hook:
	rm -rf `find $(distdir) -type d -name .svn`
	rm -rf `find $(distdir) -type d -name .git`
