/*
**  CXSC is a C++ library for eXtended Scientific Computing (V 2.5.4)
**
**  Copyright (C) 1990-2000 Institut fuer Angewandte Mathematik,
**                          Universitaet Karlsruhe, Germany
**            (C) 2000-2014 Wiss. Rechnen/Softwaretechnologie
**                          Universitaet Wuppertal, Germany   
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Library General Public
**  License as published by the Free Software Foundation; either
**  version 2 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Library General Public License for more details.
**
**  You should have received a copy of the GNU Library General Public
**  License along with this library; if not, write to the Free
**  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

##  /* CVS $Id: r_ari_x86.s,v 1.11 2014/01/30 17:23:52 cxsc Exp $ */
##  Author: Boris Hoeffgen <hoeffgen@uni-wuppertal.de>

.globl ra_addu
	.type	ra_addu, @function
ra_addu:
	pushl	%ebp
	movl	%esp, %ebp
	fstcw   -8(%ebp)             
	fstcw	-16(%ebp)
	movl	$0xB3F, -8(%ebp)
        fldcw   -8(%ebp)  
	fldl	8(%ebp)               
	faddl	16(%ebp) 
	fstpl	-24(%ebp)
	fldcw	-16(%ebp) 
	fldl	-24(%ebp)
	movl	%ebp, %esp
	popl	%ebp
	ret
	
.globl ra_addd
	.type	ra_addd, @function
ra_addd:
	pushl	%ebp
	movl	%esp, %ebp
        fstcw   -8(%ebp) 
	fstcw	-16(%ebp)       
        movl    $0x73F, -8(%ebp) 
        fldcw   -8(%ebp)        
	fldl	8(%ebp)      
	faddl	16(%ebp)
	fstpl	-24(%ebp)
	fldcw	-16(%ebp)
	fldl	-24(%ebp)
	movl	%ebp, %esp
	popl	%ebp
	ret

.globl ra_subu
	.type	ra_subu, @function
ra_subu:
	pushl	%ebp
	movl	%esp, %ebp
	fstcw   -8(%ebp)             
	fstcw	-16(%ebp)
	movl	$0xB3F, -8(%ebp)
        fldcw   -8(%ebp)  
	fldl	8(%ebp)               
	fsubl	16(%ebp) 
	fstpl	-24(%ebp)
	fldcw	-16(%ebp) 
	fldl	-24(%ebp)
	movl	%ebp, %esp
	popl	%ebp
	ret
	
.globl ra_subd
	.type	ra_subd, @function
ra_subd:
	pushl	%ebp
	movl	%esp, %ebp
        fstcw   -8(%ebp) 
	fstcw	-16(%ebp)       
        movl    $0x73F, -8(%ebp) 
        fldcw   -8(%ebp)        
	fldl	8(%ebp)      
	fsubl	16(%ebp)
	fstpl	-24(%ebp)
	fldcw	-16(%ebp)
	fldl	-24(%ebp)
	movl	%ebp, %esp
	popl	%ebp
	ret

.globl ra_muld
	.type	ra_muld, @function
ra_muld:
	pushl	%ebp
	movl	%esp, %ebp
	fstcw   -8(%ebp)             
	fstcw   -16(%ebp)             
        movl    $0x73F, -8(%ebp) 
        fldcw   -8(%ebp)        
	fldl	8(%ebp)      
	fmull	16(%ebp)
	fstpl	-24(%ebp)
	fldcw	-16(%ebp)
	fldl	-24(%ebp)
	movl	%ebp, %esp
	popl	%ebp
	ret

.globl ra_mulu
	.type	ra_mulu, @function
ra_mulu:
	pushl	%ebp
	movl	%esp, %ebp
	fstcw   -8(%ebp)             
	fstcw   -16(%ebp)             
        movl    $0xB3F, -8(%ebp) 
        fldcw   -8(%ebp)        
	fldl	8(%ebp)      
	fmull	16(%ebp)
	fstpl	-24(%ebp)
	fldcw	-16(%ebp)
	fldl	-24(%ebp)
	movl	%ebp, %esp
	popl	%ebp
	ret

.globl ra_divu
	.type	ra_divu, @function
ra_divu:
	pushl	%ebp
	movl	%esp, %ebp
	fstcw   -8(%ebp)             
	fstcw	-16(%ebp)
	movl	$0xB3F, -8(%ebp)
        fldcw   -8(%ebp)  
	fldl	8(%ebp)               
	fdivl	16(%ebp) 
	fstpl	-24(%ebp)
	fldcw	-16(%ebp) 
	fldl	-24(%ebp)
	movl	%ebp, %esp
	popl	%ebp
	ret
	
.globl ra_divd
	.type	ra_divd, @function
ra_divd:
	pushl	%ebp
	movl	%esp, %ebp
        fstcw   -8(%ebp) 
	fstcw	-16(%ebp)       
        movl    $0x73F, -8(%ebp) 
        fldcw   -8(%ebp)        
	fldl	8(%ebp)      
	fdivl	16(%ebp)
	fstpl	-24(%ebp)
	fldcw	-16(%ebp)
	fldl	-24(%ebp)
	movl	%ebp, %esp
	popl	%ebp
	ret

