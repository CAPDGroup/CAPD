#
##  fi_lib++  --- A fast interval library (Version 2.0)                     
##                                                                  
##  Copyright (C) 2001:                                                        
##                                                     
##  Werner Hofschuster, Walter Kraemer                               
##  Wissenschaftliches Rechnen/Softwaretechnologie (WRSWT)  
##  Universitaet Wuppertal, Germany                                           
##  Michael Lerch, German Tischler, Juergen Wolff von Gudenberg       
##  Institut fuer Informatik                                         
##  Universitaet Wuerzburg, Germany                                           
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
##  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
LATEX?=latex
DVIPS?=dvips
PSNUP?=psnup
DVIPDF?=dvipdf
PS2PDF14?=ps2pdf14
RM=rm -f

all:
	@${MAKE} manual2.ps
	@${MAKE} tocclean
	@${MAKE} manual2.ps
	@${MAKE} manual.pdf

manual.dvi:: manual.tex
	${LATEX} manual.tex
manual.ps:: manual.dvi
	${DVIPS} -o $@ manual.dvi
manual2.ps:: manual.ps
	${PSNUP} -2 manual.ps manual2.ps
manual.pdf:: manual.ps
	${PS2PDF14} manual.ps manual.pdf

clean:
	@${RM} manual.toc manual.aux manual.log manual.dvi manual.ps manual2.ps manual.pdf
tocclean:
	@${RM} manual.pdf manual.dvi manual.ps manual2.ps
distclean: clean tocclean
