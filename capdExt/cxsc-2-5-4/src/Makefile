#========== you shouldn't modify anything below ===========================
export FILIB = fi_lib
export RTS = rts

ifeq ($(ASM),1)
	RARI=asm/r_ari.o
	MARI=ari
endif

include rts/rts.obj
include fi_lib/fi_lib.obj
include ../CToolbox/ctoolbox.obj
include cxsc.obj

.SUFFIXES:
.SUFFIXES: .cpp .hpp .o


default:
	@echo
	@echo "C-XSC - C++ library for eXtended Scientific Computation"
	@echo "Please use the Makefile in the upper directory!"
	@echo

all: .depend libs symlinks

#test.o: test.cpp test.hpp
#	$(CXX) $(CXXFLAGS) -Wall -Winline -c test.cpp
#	$(CXX) $(CXXFLAGS) -c test.cpp

# no optimization here
test: test.o libcxsc.a
#	$(CXX) $(CXXFLAGS) -O0 -Wall -Winline test.o libcxsc.a -o test
	$(CXX) $(CXXFLAGS) test.o libcxsc.a -o test

libs: cxsc libcxsc.a libcxsc.so

ari: 
	@+cd asm; $(MAKE) ari

cxsc: 
	@+cd rts; $(MAKE) lib
	@+cd fi_lib; $(MAKE) lib
	@+cd ../CToolbox; $(MAKE) lib

libcxsc: libcxsc.a libcxsc.so

libcxsc.so: $(MARI) $(OBJ)
ifeq ($(COMPILER),ibm)
	$(CXX) $(CXXFLAGS) $(CXXOPTIMIZE) $(SHARED) $(LINKEROPTDYNLIBNAME) -o libcxsc.so.$(VER) $(OBJ_RTS) $(OBJ_CTOOLBOX) $(OBJ) $(OBJECTS) $(RARI)
else
	$(CXX) $(CXXOPTIMIZE) $(SHARED) $(LINKEROPTDYNLIBNAME) -o libcxsc.so.$(VER) $(OBJ_RTS) $(OBJ_CTOOLBOX) $(OBJ) $(OBJECTS) $(RARI)
endif

libcxsc_inl.a: $(OBJ_INL)
ifeq ($(COMPILER),sun_forte)
	$(CXX) -xar -o libcxsc_inl.a $(OBJ_INL)
else
	$(RANLIB) libcxsc_inl.a
endif

libcxsc.a: $(MARI) $(OBJ)
ifeq ($(COMPILER),sun_forte)
	$(CXX) -xar -o libcxsc.a $(OBJ_RTS) $(OBJ_CTOOLBOX) $(OBJ) $(OBJECTS) $(RARI)
else
#ifeq ($(UNAME_SYSTEM),Darwin_i386)
ifneq (,$(findstring Darwin, $(UNAME_SYSTEM))) 
	libtool -static -o libcxsc.a $(OBJ_RTS) $(OBJ_CTOOLBOX) $(OBJ) $(OBJECTS) $(RARI)
else
	$(AR) r libcxsc.a $(OBJ_RTS) $(OBJ_CTOOLBOX) $(OBJ) $(OBJECTS) $(RARI)
	$(RANLIB) libcxsc.a
endif	
endif

%.o: %.cpp
	$(CXX) -DNO_CXSC_STATUS $(CXXFLAGS) $(FPIC) $(CXXOPTIMIZE) -c $<

%.hpp: %.hpp,v
	co $@

# we had to restrict the filenames to 8+3; some say it's because
# there exists an operating system which depends on it...
symlinks:
	$(RM) cinterval.hpp intmatrix.hpp intvector.hpp l_interval.hpp
	$(RM) l_imatrix.hpp l_ivector.hpp l_rmatrix.hpp l_rvector.hpp
	$(RM) l_complex.hpp
	-$(LN) cinterva.hpp cinterval.hpp; $(LN) intmatri.hpp intmatrix.hpp
	-$(LN) intvecto.hpp intvector.hpp; $(LN) l_interv.hpp l_interval.hpp
	-$(LN) limatrix.hpp l_imatrix.hpp; $(LN) livector.hpp l_ivector.hpp
	-$(LN) lrmatrix.hpp l_rmatrix.hpp; $(LN) lrvector.hpp l_rvector.hpp
	-$(LN) l_comple.hpp l_complex.hpp;

%.h: %.h,v
	co $@

install: install_both

install_both: install_dyn install_static

install_header:
# install header files
	$(INSTDIR) $(PREFIX)/include
	$(INSTALL) ./*.hpp ./*.inl ./*.h $(PREFIX)/include
	$(INSTALL) rts/*.h $(PREFIX)/include
	$(INSTALL) fi_lib/*.hpp $(PREFIX)/include
	$(INSTALL) asm/r_ari.h $(PREFIX)/include
	$(INSTALL) ../cxscconf.h $(PREFIX)/include
	@+cd ../CToolbox; $(MAKE) install_header
	$(INSTDIR) $(PREFIX)/lib

install_dyn: install_header
# install libs
	$(RM) $(PREFIX)/liblibcxsc*
	$(INSTALL) ./libcxsc.so.$(VER) $(PREFIX)/lib
	$(LN) libcxsc.so.$(VER) libcxsc.so.2
	$(LN) libcxsc.so.2 libcxsc.so
	-cd $(PREFIX)/lib; $(LN) libcxsc.so.$(VER) libcxsc.so.2
	-cd $(PREFIX)/lib; $(LN) libcxsc.so.2 libcxsc.so

install_static: install_header
# install libs
	$(INSTALL) ./libcxsc.a $(PREFIX)/lib
#ifneq ($(COMPILER),sun_forte)
#	$(STRIP) $(PREFIX)/lib/libcxsc.a
#endif	

uninstall:
# remove libs
	$(RM) $(PREFIX)/lib/libcxsc.a $(PREFIX)/lib/libcxsc_inl.a
# remove header files
	@for i in $(UNINSTALL_OBJ); do \
	    $(RM) $(PREFIX)/include/$$i ; \
	done

clean: distclean
	@+cd rts; $(MAKE) clean
	@+cd fi_lib; $(MAKE) clean
	@+cd ../CToolbox; $(MAKE) clean
ifeq ($(COMPILER),sun_forte)
	@+cd SunWS_cache; CCadmin -clean
endif
	-$(RM) *.o libcxsc*

distclean: 
	cd .. ; \
	find . \( -name '*.orig' -o -name '*.so' -o -name '*.a' -o -name '*~' -o -name '*.bak' \
		   -o -name '#*#' -o -name '.*.orig' -o -name '*.[oa]' \
		   -o -size 0 -o -name '.*.flags' -o -name '.#*' \) \
		   -type f -print | xargs rm -f

# -M sollte jeder c++-Compiler machen, -MM spart sich die stdlib-includes
.depend depend:
ifeq ($(COMPILER),sun_forte)
	umask 002; $(CXX) $(CXXFLAGS) $(DEPENDOPT) $(SRCS) $(VSRCS) > ".depend"
else
	umask 002; $(CXX) $(DEPENDOPT) $(CXXFLAGS) $(SRCS) $(VSRCS) > ".depend"
endif

.PHONY: symlinks clean depend all libs install uninstall

include .depend

#
#	create spec-file for rpm 
#
spec:
	. ../scripts/mkspec >../cxsc.spec

#
#	Build a tar ball, generate an rpm from it and pack the result
#
rpm:	clean spec
	set -e ; \
	cd $(TOPDIR)/.. ; \
	ln -sf $(TOPDIR) $(CXSCPATH) ; \
	tar -cvz --exclude CVS -f $(CXSCPATH).tar.gz $(CXSCPATH)/. ; \
	$(RM) $(CXSCPATH) ; \
	cd $(TOPDIR) ; \
	rpm -ta $(TOPDIR)/../$(CXSCPATH).tar.gz ; \
	mv /usr/src/packages/RPMS/i386/$(CXSCPATH)-$(RELEASE).i386.rpm $(TOPDIR)/$(CXSCPATH)-$(RELEASE).rpm ; \
	$(RM) $(TOPDIR)/../$(CXSCPATH).tar.gz ; \
	$(RM) $(TOPDIR)/cxsc.spec ; \
	$(RM) $(CXSCPATH) ; \
	$(RM) /usr/src/packages/SRPMS/$(CXSCPATH)-$(RELEASE).src.rpm ; \
	$(RMDIR) /usr/src/packages/BUILD/$(CXSCPATH) ; \
	$(RMDIR) /usr/src/packages/SOURES/$(CXSCPATH) ; \
	$(RM) /usr/src/packages/SPECS/cxsc.spec ; \

#
#	Build and generate a sunpkg 
#
sunpkg: clean all
	$(INSTDIR) $(TMPDIR)/lib ; \
	$(INSTALL) ./libcxsc.a $(TMPDIR)/lib/ ; \
	$(INSTALL) ./libcxsc_inl.a $(TMPDIR)/lib/ ; \
	$(INSTDIR) $(TMPDIR)/include ; \
	$(INSTALL) ./*.hpp ./*.inl ./*.h $(TMPDIR)/include ; \
	$(INSTALL) rts/*.h $(TMPDIR)/include ; \
	$(INSTALL) fi_lib/fi_lib.hpp $(TMPDIR)/include ; \
	$(INSTALL) ../cxscconf.h $(TMPDIR)/include ; \
	cd $(TMPDIR) ; \
	find . -print | pkgproto >out ; \
	awk -f ../scripts/prototype.awk out >prototype ; \
	. ../scripts/mkpkginfo >$(TMPDIR)/pkginfo ; \
	pkgmk -r `pwd` ; \
	cd $(BUILDTMP) ; \
	pkgtrans -s `pwd` /tmp/$(CXSCPATH)-$(RELEASE).$(ARCH).pkg $(PKGNAME) ; \
	$(RMDIR) $(TMPDIR) ; \
	$(RMDIR) $(BUILDTMP)/$(PKGNAME) ; \
	gzip /tmp/$(CXSCPATH)-$(RELEASE).$(ARCH).pkg ; \

#
#	generate a backup tar archive
#
backup:	clean
	set -e ; \
	cd $(TOPDIR)/.. ; \
	ln -sf $(TOPDIR) $(CXSCPATH) ; \
	tar -cvf $(CXSCPATH)-$(RELEASE).tar $(CXSCPATH)/. ; \
	gzip $(CXSCPATH)-$(RELEASE).tar ; \
	rm $(CXSCPATH) ; \
