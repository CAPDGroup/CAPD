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

##  /* CVS $Id: r_ari_x86_i.s,v 1.11 2014/01/30 17:23:52 cxsc Exp $ */
##  Author: Boris Hoeffgen <hoeffgen@uni-wuppertal.de>

.globl ra_addu
	.type	ra_addu, @function
ra_addu:
	pushl	%ebp
	movl	%esp, %ebp
	stmxcsr -8(%ebp)             
	stmxcsr -16(%ebp)             
	movl	$0x5F80, -8(%ebp)
        ldmxcsr -8(%ebp)  
	addsd	%xmm1, %xmm0
	ldmxcsr -16(%ebp)
	movl	%ebp, %esp
	popl	%ebp
	ret
	
.globl ra_addd
	.type	ra_addd, @function
ra_addd:
	pushl	%ebp
	movl	%esp, %ebp
	stmxcsr -8(%ebp)             
	stmxcsr -16(%ebp)             
	movl	$0x3F80, -8(%ebp)
        ldmxcsr -8(%ebp)  
    	addsd   %xmm1,%xmm0 
	ldmxcsr -16(%ebp)
	movl	%ebp, %esp
	popl	%ebp
	ret

.globl ra_subu
	.type	ra_subu, @function
ra_subu:
	pushl	%ebp
	movl	%esp, %ebp
	stmxcsr -8(%ebp)             
	stmxcsr -16(%ebp)             
	movl	$0x5F80, -8(%ebp)
        ldmxcsr -8(%ebp)  
	subsd	%xmm1, %xmm0
	ldmxcsr -16(%ebp)
	movl	%ebp, %esp
	popl	%ebp
	ret
	
.globl ra_subd
	.type	ra_subd, @function
ra_subd:
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

.globl ra_muld
	.type	ra_muld, @function
ra_muld:
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

.globl ra_mulu
	.type	ra_mulu, @function
ra_mulu:
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

.globl ra_divu
	.type	ra_divu, @function
ra_divu:
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
	
.globl ra_divd
	.type	ra_divd, @function
ra_divd:
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

