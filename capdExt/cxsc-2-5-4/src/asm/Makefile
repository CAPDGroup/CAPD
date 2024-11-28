##
##  CXSC is a C++ library for eXtended Scientific Computing
##
##  Copyright (C) 1990-2000 Institut fuer Angewandte Mathematik,
##                          Universitaet Karlsruhe, Germany
##            (C) 2000-2005 Wiss. Rechnen/Softwaretechnologie
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

ARCH?=$(shell uname -m | sed -e s/i.86/i386/ -e s/sun4u/sparc/)

ifeq "$(ARCH)" "x86_64"
	ifeq "$(BIT)" "32"
		ARI=r_ari_x86
	else
		ARI=r_ari_x86_64
	endif
endif

ifeq "$(ARCH)" "i386"
	ARI=r_ari_x86
endif

ifeq ($(UNAME_SYSTEM),Darwin_i386)
	ifeq "$(BIT)" "32"
		ARI=r_ari_x86_mac
	else
		ARI=r_ari_x86_64_mac
	endif
endif

.SUFFIXES: .s .o

ari: $(ARI)
	@mv $(ARI).o r_ari.o	

$(ARI): $(ARI).s
	$(CC) -c $(CCFLAGS) $(FPIC) $(CCOPTIMIZE) $<

clean:
	rm -f *.o
	
