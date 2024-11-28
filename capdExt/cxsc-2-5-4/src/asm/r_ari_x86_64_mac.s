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

##  /* CVS $Id: r_ari_x86_64_mac.s,v 1.10 2014/01/30 17:23:52 cxsc Exp $ */
##  Author: Boris Hoeffgen <hoeffgen@uni-wuppertal.de>

.globl _ra_addu
_ra_addu:
	pushq	%rbp
	movq	%rsp, %rbp
	stmxcsr -8(%rbp)             
	stmxcsr -16(%rbp)             
	movl	$0x5F80, -8(%rbp)
        ldmxcsr -8(%rbp)  
	addsd	%xmm1, %xmm0
	ldmxcsr -16(%rbp)
	movq	%rbp, %rsp
	popq	%rbp
	ret
	
.globl _ra_addd
_ra_addd:
	pushq	%rbp
	movq	%rsp, %rbp
	stmxcsr -8(%rbp)             
	stmxcsr -16(%rbp)             
	movl	$0x3F80, -8(%rbp)
        ldmxcsr -8(%rbp)  
    	addsd   %xmm1,%xmm0 
	ldmxcsr -16(%rbp)
	movq	%rbp, %rsp
	popq	%rbp
	ret

.globl _ra_subu
_ra_subu:
	pushq	%rbp
	movq	%rsp, %rbp
	stmxcsr -8(%rbp)             
	stmxcsr -16(%rbp)             
	movl	$0x5F80, -8(%rbp)
        ldmxcsr -8(%rbp)  
	subsd	%xmm1, %xmm0
	ldmxcsr -16(%rbp)
	movq	%rbp, %rsp
	popq	%rbp
	ret
	
.globl _ra_subd
_ra_subd:
	pushq	%rbp
	movq	%rsp, %rbp
	stmxcsr -8(%rbp)             
	stmxcsr -16(%rbp)             
	movl	$0x3F80, -8(%rbp)
        ldmxcsr -8(%rbp)  
    	subsd   %xmm1,%xmm0 
	ldmxcsr -16(%rbp)
	movq	%rbp, %rsp
	popq	%rbp
	ret

.globl _ra_muld
_ra_muld:
	pushq	%rbp
	movq	%rsp, %rbp
	stmxcsr -8(%rbp)             
	stmxcsr -16(%rbp)             
	movl	$0x3F80, -8(%rbp)
        ldmxcsr -8(%rbp)  
    	mulsd   %xmm1,%xmm0 
	ldmxcsr -16(%rbp)
	movq	%rbp, %rsp
	popq	%rbp
	ret

.globl _ra_mulu
_ra_mulu:
	pushq	%rbp
	movq	%rsp, %rbp
	stmxcsr -8(%rbp)             
	stmxcsr -16(%rbp)             
	movl	$0x5F80, -8(%rbp)
        ldmxcsr -8(%rbp)  
	mulsd	%xmm1, %xmm0
	ldmxcsr -16(%rbp)
	movq	%rbp, %rsp
	popq	%rbp
	ret

.globl _ra_divu
_ra_divu:
	pushq	%rbp
	movq	%rsp, %rbp
	stmxcsr -8(%rbp)             
	stmxcsr -16(%rbp)             
	movl	$0x5F80, -8(%rbp)
        ldmxcsr -8(%rbp)  
	divsd	%xmm1, %xmm0
	ldmxcsr -16(%rbp)
	movq	%rbp, %rsp
	popq	%rbp
	ret
	
.globl _ra_divd
_ra_divd:
	pushq	%rbp
	movq	%rsp, %rbp
	stmxcsr -8(%rbp)             
	stmxcsr -16(%rbp)             
	movl	$0x3F80, -8(%rbp)
        ldmxcsr -8(%rbp)  
    	divsd   %xmm1,%xmm0 
	ldmxcsr -16(%rbp)
	movq	%rbp, %rsp
	popq	%rbp
	ret

