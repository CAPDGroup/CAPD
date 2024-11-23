#
# This file is the Makefile for the C-XSC Toolbox
# ===============================================
#
##
##  CXSC is a C++ library for eXtended Scientific Computing 
##
##  Copyright (C) 1990-2000 Institut fuer Angewandte Mathematik,
##                          Universitaet Karlsruhe, Germany
##            (C) 2000-2009 Wiss. Rechnen/Softwaretechnologie
##                          Universitaet Wuppertal, Germany
##
##  This library is free software; you can redistribute it and/or
##  modify it under the terms of the GNU Library General Public
##  License as published by the Free Software Foundation; either
##  version 2 of the License, or (at your option) any later version.
##
##  This library is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
##  Library General Public License for more details.
##
##  You should have received a copy of the GNU Library General Public
##  License along with this library; if not, write to the Free
##  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
##

#
#=================== don't modify the following values ========================#
#
include ctoolbox.obj
#the C-XSC root directory
CXSCDIR = $(PREFIX)
export CXSCDIR

# Settings for all systems ---------------------------------------------------
M = Modules#               		the toolbox modules directory
P = Programs#              		the toolbox programs directory

CXXINC=-I$(CXSCDIR)/include -I$(M) -I../src -I../src/rts -I../src/fi_lib -I../src/asm -I../ -I./#     additional include path
CXXFLAGS=$(CXXINC) $(CXXOPTS)#
ifeq ($(DYN),1)
ifeq ($(UNAME_SYSTEM),Darwin_i386)
	RPATH = $(CXSCDIR)/lib/libcxsc.so
else
	RPATH = -Wl,-R$(CXSCDIR)/lib
	export RPATH
endif
endif
LIBRARIES=-L$(CXSCDIR)/lib $(LIBS)

# Suffix rules
#-------------
.SUFFIXES:
.SUFFIXES: .cpp .o .a 

.cpp.o:
	$(CXX) -DNO_CXSC_STATUS -c $(FPIC) $(CXXFLAGS) $(CXXOPTIMIZE) $*.cpp -o $*.o

.cpp.a:
	$(CXX) -DNO_CXSC_STATUS -c $(FPIC) $(CXXFLAGS) $(CXXOPTIMIZE) $*.cpp -o $*.o
	@rm -f $*.o


# Default rule
#-------------

default:
	@echo
	@echo "C++ Toolbox for Verified Computing"
	@echo "you additionally need the C-XSC library"
	@echo

# Rules for generating executables of the sample programs
#--------------------------------------------------------
PRJ1  = cpz_ex
PRJ2  = ddf_ex
PRJ3  = gop_ex
PRJ4  = gop1_ex
PRJ5  = hess_ex
PRJ6  = jac_ex
PRJ7  = lop_ex
PRJ8  = lss_ex
PRJ9  = nlfz_ex
PRJ10 = nlss_ex
PRJ11 = rpe_ex
PRJ12 = xev_ex1
PRJ13 = xev_ex2
PRJ14 = xev_ex3
PRJ15 = fastlss_ex

$(PRJ1): $(P)/$(PRJ1).cpp
	$(CXX) -o $(P)/$@ $(P)/$@.cpp $(CXXFLAGS) $(CXXOPTIMIZE) $(RPATH) $(LIBRARIES)

$(PRJ2): $(P)/$(PRJ2).cpp
	$(CXX) -o $(P)/$@ $(P)/$@.cpp $(CXXFLAGS) $(CXXOPTIMIZE) $(RPATH) $(LIBRARIES)

$(PRJ3): $(P)/$(PRJ3).cpp
	$(CXX) -o $(P)/$@ $(P)/$@.cpp $(CXXFLAGS) $(CXXOPTIMIZE) $(RPATH) $(LIBRARIES)

$(PRJ4): $(P)/$(PRJ4).cpp
	$(CXX) -o $(P)/$@ $(P)/$@.cpp $(CXXFLAGS) $(CXXOPTIMIZE) $(RPATH) $(LIBRARIES)

$(PRJ5): $(P)/$(PRJ5).cpp
	$(CXX) -o $(P)/$@ $(P)/$@.cpp $(CXXFLAGS) $(CXXOPTIMIZE) $(RPATH) $(LIBRARIES)

$(PRJ6): $(P)/$(PRJ6).cpp
	$(CXX) -o $(P)/$@ $(P)/$@.cpp $(CXXFLAGS) $(CXXOPTIMIZE) $(RPATH) $(LIBRARIES)

$(PRJ7): $(P)/$(PRJ7).cpp
	$(CXX) -o $(P)/$@ $(P)/$@.cpp $(CXXFLAGS) $(CXXOPTIMIZE) $(RPATH) $(LIBRARIES)
 
$(PRJ8): $(P)/$(PRJ8).cpp
	$(CXX) -o $(P)/$@ $(P)/$@.cpp $(CXXFLAGS) $(CXXOPTIMIZE) $(RPATH) $(LIBRARIES)

$(PRJ9): $(P)/$(PRJ9).cpp
	$(CXX) -o $(P)/$@ $(P)/$@.cpp $(CXXFLAGS) $(CXXOPTIMIZE) $(RPATH) $(LIBRARIES)

$(PRJ10): $(P)/$(PRJ10).cpp
	$(CXX) -o $(P)/$@ $(P)/$@.cpp $(CXXFLAGS) $(CXXOPTIMIZE) $(RPATH) $(LIBRARIES)

$(PRJ11): $(P)/$(PRJ11).cpp
	$(CXX) -o $(P)/$@ $(P)/$@.cpp $(CXXFLAGS) $(CXXOPTIMIZE) $(RPATH) $(LIBRARIES)

$(PRJ12): $(P)/$(PRJ12).cpp
	$(CXX) -o $(P)/$@ $(P)/$@.cpp $(CXXFLAGS) $(CXXOPTIMIZE) $(RPATH) $(LIBRARIES)

$(PRJ13): $(P)/$(PRJ13).cpp
	$(CXX) -o $(P)/$@ $(P)/$@.cpp $(CXXFLAGS) $(CXXOPTIMIZE) $(RPATH) $(LIBRARIES)

$(PRJ14): $(P)/$(PRJ14).cpp
	$(CXX) -o $(P)/$@ $(P)/$@.cpp $(CXXFLAGS) $(CXXOPTIMIZE) $(RPATH) $(LIBRARIES)

$(PRJ15): $(P)/$(PRJ15).cpp
	$(CXX) -o $(P)/$@ $(P)/$@.cpp $(CXXFLAGS) $(CXXOPTIMIZE) $(RPATH) $(LIBRARIES)

ALL_PROGRAMS = $(PRJ1)  $(PRJ2)  $(PRJ3)  $(PRJ4)  $(PRJ5)  $(PRJ6)  \
               $(PRJ7)  $(PRJ8)  $(PRJ9)  $(PRJ10) $(PRJ11) $(PRJ12) \
               $(PRJ13) $(PRJ14) $(PRJ15)

all: $(ALL_PROGRAMS)	
	@echo "all toolbox sample programs are (now) up-to-date"
	@echo

# Rules for generating the test programs
#---------------------------------------
TDIR = test
TPROG1  = tester
TPROG2  = sumup

$(TPROG1): $(TDIR)/$(TPROG1).cpp
	$(CXX) -o $(TDIR)/$(TPROG1) $(CXXINC) $(CXXOPTS) $(CXXOPTIMIZE) $(TDIR)/$(TPROG1).cpp $(RPATH) $(LIBRARIES)

$(TPROG2): $(TDIR)/$(TPROG2).cpp
	$(CXX) -I../ -o $(TDIR)/$(TPROG2) $(TDIR)/$(TPROG2).cpp 

lib: $(OBJ_CTOOLBOX)
	

test: $(TPROG1) $(TPROG2)
	($(P)/cpz_ex < $(P)/cpz_ex.in1) > $(TDIR)/cpz_ex.out1
	($(P)/cpz_ex < $(P)/cpz_ex.in2) > $(TDIR)/cpz_ex.out2
	($(P)/gop1_ex < $(P)/gop1_ex.in) > $(TDIR)/gop1_ex.out
	($(P)/gop_ex < $(P)/gop_ex.in) > $(TDIR)/gop_ex.out
	($(P)/hess_ex < $(P)/hess_ex.in) > $(TDIR)/hess_ex.out
	($(P)/jac_ex < $(P)/jac_ex.in) > $(TDIR)/jac_ex.out
	($(P)/lop_ex < $(P)/lop_ex.in1) > $(TDIR)/lop_ex.out1
	($(P)/lop_ex < $(P)/lop_ex.in2) > $(TDIR)/lop_ex.out2
	($(P)/lss_ex < $(P)/lss_ex.in1) > $(TDIR)/lss_ex.out1
	($(P)/lss_ex < $(P)/lss_ex.in2) > $(TDIR)/lss_ex.out2
	($(P)/nlfz_ex < $(P)/nlfz_ex.in) > $(TDIR)/nlfz_ex.out
	($(P)/nlss_ex < $(P)/nlss_ex.in1) > $(TDIR)/nlss_ex.out1
	($(P)/nlss_ex < $(P)/nlss_ex.in2) > $(TDIR)/nlss_ex.out2
	($(P)/rpe_ex < $(P)/rpe_ex.in1) > $(TDIR)/rpe_ex.out1
	($(P)/rpe_ex < $(P)/rpe_ex.in2) > $(TDIR)/rpe_ex.out2
	($(P)/xev_ex1 < $(P)/xev_ex1.in) > $(TDIR)/xev_ex1.out
	($(P)/xev_ex2 < $(P)/xev_ex2.in) > $(TDIR)/xev_ex2.out
	($(P)/xev_ex3 < $(P)/xev_ex3.in) > $(TDIR)/xev_ex3.out
	rm -f $(TDIR)/ergebnis.dat
	touch $(TDIR)/ergebnis.dat
	cd $(TDIR); ./tester cpz_ex.out1 cpz_ex.vgl1
	cd $(TDIR); ./tester cpz_ex.out2 cpz_ex.vgl2
	cd $(TDIR); ./tester gop1_ex.out gop1_ex.vgl
	cd $(TDIR); ./tester gop_ex.out gop_ex.vgl
	cd $(TDIR); ./tester hess_ex.out hess_ex.vgl
	cd $(TDIR); ./tester jac_ex.out jac_ex.vgl
	cd $(TDIR); ./tester lop_ex.out1 lop_ex.vgl1
	cd $(TDIR); ./tester lop_ex.out2 lop_ex.vgl2
	cd $(TDIR); ./tester lss_ex.out1 lss_ex.vgl1
	cd $(TDIR); ./tester lss_ex.out2 lss_ex.vgl2
	cd $(TDIR); ./tester nlfz_ex.out nlfz_ex.vgl
	cd $(TDIR); ./tester nlss_ex.out1 nlss_ex.vgl1
	cd $(TDIR); ./tester nlss_ex.out2 nlss_ex.vgl2
	cd $(TDIR); ./tester rpe_ex.out1 rpe_ex.vgl1
	cd $(TDIR); ./tester rpe_ex.out2 rpe_ex.vgl2
	cd $(TDIR); ./tester xev_ex1.out xev_ex1.vgl
	cd $(TDIR); ./tester xev_ex2.out xev_ex2.vgl
	cd $(TDIR); ./tester xev_ex3.out xev_ex3.vgl
	cd $(TDIR); ./sumup

install_header:
# install header files
	$(INSTDIR) $(PREFIX)/include
	$(INSTALL) $(M)/*.hpp  $(PREFIX)/include
	$(INSTALL) $(M)/*.inl  $(PREFIX)/include

install: install_header
	$(INSTDIR) $(PREFIX)/examples
	@cd Programs; $(INSTALL) cpz_ex ddf_ex gop_ex gop1_ex hess_ex jac_ex lop_ex lss_ex nlfz_ex nlss_ex rpe_ex xev_ex1 xev_ex2 xev_ex3 fastlss_ex $(PREFIX)/examples
	@cd Programs; $(INSTALL) cpz_ex.cpp ddf_ex.cpp gop_ex.cpp gop1_ex.cpp hess_ex.cpp jac_ex.cpp lop_ex.cpp lss_ex.cpp nlfz_ex.cpp nlss_ex.cpp rpe_ex.cpp xev_ex1.cpp xev_ex2.cpp xev_ex3.cpp fastlss_ex.cpp $(PREFIX)/examples

installprog:
	$(INSTDIR) $(PREFIX)/examples
	@cd Programs; $(INSTALLEXE) cpz_ex ddf_ex gop_ex gop1_ex hess_ex jac_ex lop_ex lss_ex nlfz_ex nlss_ex rpe_ex xev_ex1 xev_ex2 xev_ex3 fastlss_ex $(PREFIX)/examples
	@cd Programs; $(INSTALL) cpz_ex.cpp ddf_ex.cpp gop_ex.cpp gop1_ex.cpp hess_ex.cpp jac_ex.cpp lop_ex.cpp lss_ex.cpp nlfz_ex.cpp nlss_ex.cpp rpe_ex.cpp xev_ex1.cpp xev_ex2.cpp xev_ex3.cpp fastlss_ex.cpp $(PREFIX)/examples

uninstall:
# remove header files
	@for i in $(UNINSTALL_OBJ); do \
		$(RM) $(PREFIX)/include/$$i ; \
	done

clean:
	rm -f $(P)/$(PRJ1) $(P)/$(PRJ2) $(P)/$(PRJ3) $(P)/$(PRJ4) \
	$(P)/$(PRJ5) $(P)/$(PRJ6) $(P)/$(PRJ7) $(P)/$(PRJ8) $(P)/$(PRJ9) \
	$(P)/$(PRJ10) $(P)/$(PRJ11) $(P)/$(PRJ12) $(P)/$(PRJ13) \
	$(P)/$(PRJ14) $(P)/$(PRJ15)
	rm -f $(P)/*.o
	rm -f test/sumup test/tester test/*.out*

.phony: test
