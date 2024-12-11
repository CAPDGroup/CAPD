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

/* CVS $Id: r_ftrp.c,v 1.21 2014/01/30 17:24:12 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : r_ftrp.c                              */
/*                                                              */
/*      Entries         : void r_ftrp(...)                      */
/*                        ... arguments depends on compilersw.  */
/*                                                              */
/*      Description     : floating point trap handler for       */
/*                        Pascal-XSC arithmetic.                */
/*                                                              */
/*      Note            : Version for SUN4_OS4_C                */
/*                                                              */
/*                        r_ftrp =b=                            */
/*                        add dummy entry for non IEEEE_HW      */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
#endif

#ifdef IEEE_HARDWARE
#if SUN4_OS4_C+SUN4_GNU_C+SUN4_CPP_C

/* this have already been done in r_fcth.h =b=
 * #include <signal.h>
 * #include <sys/ieeefp.h>
 * #include <floatingpoint.h>
 */

#ifdef LINT_ARGS
void r_ftrp( int sig, int code, struct sigcontext *scp , char *addr)
#else
void r_ftrp(sig, code, scp, addr) 
int sig;
int code;
struct sigcontext *scp ;
char *addr;
#endif
{
	fprintf(stderr, "ieee exception code %x occured at pc %x\n",
		code, scp->sc_pc);

        /* test status flag for exception occurrence */
        switch (code)
	{
     
        /* test status flag for overflow occurrence */
           case FPE_FLTOVF_TRAP : 
                       if (e_of_e()) 
                          e_trap(OVERFLOW,0);
		       else
                          e_sofo();
                       break ;

        /* test status flag for underflow occurrence */
           case FPE_FLTUND_TRAP : 
                       if (e_uf_e()) 
                          e_trap(UNDERFLOW,0);
		       else
		          e_sufo();
                       break ;

        /* test status flag for invalid-operation occurrence */
           case FPE_FLTOPERR_TRAP :
                       if (e_io_e()) 
                          e_trap(INV_OP, 2, E_TMSG, 69);
		       else
                          e_sioo();
                       break ;

        /* test status flag for division-by-zero occurrence */
           case FPE_FLTDIV_TRAP :
                       if (e_dz_e()) 
                          e_trap(DIV_BY_ZERO,0);
		       else
                          e_sdzo();
                       break ;

        /* test status flag for division-by-zero occurrence */
           case FPE_FLTINEX_TRAP :
                       if (e_ie_e()) 
                          e_trap(INEXACT,0);
		       else
                          e_sieo();
                       break ;
       
           }
}

#endif /* SUN4_OS4_C+SUN4_GNU_C+SUN4_CPP_C */
#if SUN4_OS5_GNU_C

/* current exceptions */
#define EXC_NVA         0x00000010      /* invalid operation */
#define EXC_OFA         0x00000008      /* overflow */
#define EXC_UFA         0x00000004      /* underflow */
#define EXC_DZA         0x00000002      /* underflow */
#define EXC_NXA         0x00000001      /* inexact */

#ifdef LINT_ARGS
void r_ftrp(int sig, siginfo_t *sip, ucontext_t *ucp)
#else
void r_ftrp(sig, sip, ucp) 
int sig;
siginfo_t *sip;
ucontext_t *ucp;
#endif
{
        int code, op, fsr, fsr1;
	double op1, op2, result;
	
	fsr = ucp->uc_mcontext.fpregs.fpu_fsr;
	code = fsr & 0x01f;
	fprintf(stderr, "ieee exception code %x occured at pc %x\n",
		code, ucp->uc_mcontext.gregs[REG_PC]);

	/* get fpu op code */
	if (ucp->uc_mcontext.fpregs.fpu_qcnt != 1) {
	    fprintf(stderr, "cannot determine floating-point op code\n");
	    e_trap(INV_OP+E_EXIT, 0);
	}
	op = ucp->uc_mcontext.fpregs.fpu_q->FQu.fpq.fpq_instr;
	/* 81a08844 add, 81a088c4 sub, 81a08944 mul, 81a089c4 div */

	/* get operands and determine default result */
	op1 = ucp->uc_mcontext.fpregs.fpu_fr.fpu_dregs[1];
	op2 = ucp->uc_mcontext.fpregs.fpu_fr.fpu_dregs[2];
	fsr1 = fsr & 0xf0700fff; 	
	r_lfsx(fsr1);			/* disable all traps */
	switch (op)
	{
	   case 0x81a08844 : result = op1+op2; break;
	   case 0x81a088c4 : result = op1-op2; break;
	   case 0x81a08944 : result = op1*op2; break;
	   case 0x81a089c4 : result = op1/op2; break;
	   default :	fprintf(stderr, "unknown floating-point op code\n");
	   		e_trap(INV_OP+E_EXIT, 0);
			break;
	}
	r_lfsx(fsr);			/* restore fsr */

        /* test status flag for exception occurrence */
	/* set pxsc flags or generate pxsc exception */
	if (code & EXC_NVA) {
	    if (e_io_e())
		e_trap(INV_OP+E_IEEE, 8, E_TMSG, 69,
			E_TDBL+E_TEXT(1), &op1, E_TDBL+E_TEXT(2), &op2,
			E_TDBL+E_TEXT(3), &result);
	    else
		e_sioo();
	}
	
	if (code & EXC_OFA) {
	    if (e_of_e())
		e_trap(OVERFLOW+E_IEEE, 6,
			E_TDBL+E_TEXT(1), &op1, E_TDBL+E_TEXT(2), &op2,
			E_TDBL+E_TEXT(3), &result);
	    else
		e_sofo();
	}
	
	if (code & EXC_UFA) {
	    if (e_uf_e())
		e_trap(UNDERFLOW+E_IEEE, 6,
			E_TDBL+E_TEXT(1), &op1, E_TDBL+E_TEXT(2), &op2,
			E_TDBL+E_TEXT(3), &result);
	    else
		e_sufo();
	}
	
	if (code & EXC_DZA) {
	    if (e_dz_e())
		e_trap(DIV_BY_ZERO+E_IEEE, 6,
			E_TDBL+E_TEXT(1), &op1, E_TDBL+E_TEXT(2), &op2,
			E_TDBL+E_TEXT(3), &result);
	    else
		e_sdzo();
	}
	
	if (code & EXC_NXA) {
	    if (e_ie_e())
		e_trap(INEXACT, 0);
	    else
		e_sieo();
	}

	/* store result */
	ucp->uc_mcontext.fpregs.fpu_fr.fpu_dregs[0] = result;

        signal(SIGFPE, (void (*)(int)) r_ftrp);	 /* restore signal handler */
}

#endif /* SUN4_OS5_GNU_C */
#if HP_9000_C

/* exceptions */
#define EXC_NVA		0x80000000	/* invalid operation */
#define EXC_DZA		0x40000000	/* div. by zero */
#define EXC_OFA		0x20000000	/* overflow */
#define EXC_UFA		0x10000000	/* underflow */
#define EXC_NXA		0x08000000	/* inexact */

#ifdef LINT_ARGS
void r_ftrp(int code)
#else
void r_ftrp(code)
int code;
#endif
{
	fprintf(stderr, "ieee exception code %x occured\n", code);

	/* set pxsc flags or generate pxsc exception */
	if (code & EXC_NVA)
	    if (e_io_e())
		e_trap(INV_OP, 2, E_TMSG, 69);
	    else
		e_sioo();
	
	else if (code & EXC_OFA)
	    if (e_of_e())
		e_trap(OVERFLOW, 0);
	    else
		e_sofo();
	
	else if (code & EXC_UFA)
	    if (e_uf_e())
		e_trap(UNDERFLOW, 0);
	    else
		e_sufo();
	
	else if (code & EXC_DZA)
	    if (e_dz_e())
		e_trap(DIV_BY_ZERO, 0);
	    else
		e_sdzo();
	
	else if (code & EXC_NXA)
	    if (e_ie_e())
		e_trap(INEXACT, 0);
	    else
		e_sieo();
}

#endif /* HP_9000_C */
#if IBM_LINUX_C+IBM_EMX_C

/* Exceptions koennen folgendermassen auftreten:                      */
/* EXC_NVA (invalid)   bei jeder Operation mit qNaNs, ausserdem bei   */
/*                     inf-inf, 0*inf, inf/inf, 0/0 und bei FP        */
/*                     stack over/underflow                           */
/* EXC_OFA (overflow)  nur bei fst (store)                            */
/* EXC_UFA (underflow) nur bei fst (store)                            */
/* EXC_DZA (zero div)  bei der Division                               */
/* EXC_NXA (inexact)   bei +,-,*,/,fst                                */

/* Behandlung, falls das Programm fortgesetzt wird:                   */
/* EXC_NXA:  ein FP stack over/underflow beendet das Programm.        */
/*           ansonsten wird bei +,-,*,/ eine qNaN in %st(0) zurueck-  */
/*           gegeben, bei fst eine qNaN an der Zieladresse abgelegt   */
/* EXC_OFA:  der entsprechende Rueckgabewert wird an der Zieladresse  */
/*           abgelegt                                                 */
/* EXC_UFA:  wie EXC_OFA                                              */
/* EXC_DZA:  ein entsprechendes Unendlich wird in %st(0) abgelegt     */
/* EXC_NXA:  bei +,-,*,/ wird das Ergebnis wieder in %st(0) gebracht, */
/*           bei fst wird gerundet und ans Ziel gespeichert           */

/* exception flags */
#define EXC_NVA         0x0001          /* invalid operation */
#define EXC_OFA         0x0008          /* overflow */
#define EXC_UFA         0x0010          /* underflow */
#define EXC_DZA         0x0004          /* underflow */
#define EXC_NXA         0x0020          /* inexact */
#define EXC_STK		0x0040		/* stack overflow/underflow */

/* coprocessor state */
struct i387_hard_struct {
        long    cwd;		/* control word */
        long    swd;		/* status ord - cleared by trap handler */
        long    twd;		/* tag word - cleared by trap handler */
        long    fip;		/* instruction pointer */
        long    fcs;		/* op code (hi) and orig. status word (low) */
        long    foo;		/* operand offset */
        long    fos;		/* orig. tag word */
        long    st[20];         /* 8*10 bytes for each FP-reg = 80 bytes */
};

static struct i387_hard_struct fpu;

#ifdef LINT_ARGS
void r_ftrp(int sig)
#else
void r_ftrp(sig)
int sig;
#endif
{
	a_real op1, op2, result, *resp;
	unsigned short status, control, code;
	int loadres = 0;	/* load result onto fp stack ? */
	int storeres = 0;	/* store result at destination address */
	int type = 0;		/* 1=load, 2=store, 3=operation */

        __asm__("fnsave %0" : "=m" (fpu));
        status = fpu.fcs & 0x0000ffff;
	control = (fpu.cwd & 0x0f00) | 0x03f;	/* all exceptions masked */
	__asm__("fldcw %0" : "=m" (control));	/* during exc. routine */

	fprintf(stderr, "ieee exception code %x,%x occured at pc %x\n",
		status, fpu.cwd&0x0ffff, fpu.fip);

	/* handle stack overflow/underflow -> exit program */
	if (status & EXC_STK) {
	    fprintf(stderr, "floating point stack overflow/underflow\n");
	    e_trap(INV_OP+E_EXIT, 0);
	}

        /* determine if this was an operation or a store */
	code = (fpu.fcs >> 16) & 0x07ff;	/* get fpu op code */
	switch (code & 0x0738)			/* zero out mod and r/m */
	{   case 0x500:	type=1; break;	/* fldl */
	    case 0x510:			/* fstl */
	    case 0x518:	type=2; break;	/* fstpl */
	    case 0x400:			/* fadd */
	    case 0x420:			/* fsub */
	    case 0x408:			/* fmul */
	    case 0x430:	type=3; break;	/* fdiv */
	}

	/* handle stack overflow/underflow -> exit program */
	if (!type) {
	    fprintf(stderr, "unknown fp op code in exception handler\n");
	    e_trap(INV_OP+E_EXIT, 0);
	}

#ifdef IEEE_DEBUG
	fprintf(stderr, "ieee exception opcode is %x\n", code&0x738);
	fprintf(stderr, "handler: sig %x cwd %x swd %x twd %x\n", sig,fpu.cwd,fpu.swd,fpu.twd);
	fprintf(stderr, "handler: fip %x fcs %x foo %x fos %x\n", fpu.fip,fpu.fcs,fpu.foo,fpu.fos);
	fprintf(stderr, "handler: st 1 %08.8x %08.8x %08.8x %08.8x %08.8x\n", fpu.st[0], fpu.st[1], fpu.st[2], fpu.st[3], fpu.st[4]);
	fprintf(stderr, "handler: st 2 %08.8x %08.8x %08.8x %08.8x %08.8x\n", fpu.st[5], fpu.st[6], fpu.st[7], fpu.st[8], fpu.st[9]);
	fprintf(stderr, "handler: st 3 %08.8x %08.8x %08.8x %08.8x %08.8x\n", fpu.st[10], fpu.st[11], fpu.st[12], fpu.st[13], fpu.st[14]);
	fprintf(stderr, "handler: st 4 %08.8x %08.8x %08.8x %08.8x %08.8x\n", fpu.st[15], fpu.st[16], fpu.st[17], fpu.st[18], fpu.st[19]);
	fprintf(stderr, "handler: *oos %08.8x %08.8x\n\n", *((long *)fpu.foo), *(((long *)fpu.foo)+1));
#endif /* IEEE_DEBUG */

	/* get original operands */
	switch (type)
	{   case 1:	op1 = (a_real)(*((double *)fpu.foo));	/* fldl */
			loadres = 1;
			break;
	    case 2:	__asm__("fldt %0" :  : "m" (fpu.st[0])); /* fstl */
	    		__asm__("fstpl %0" : "=m" (op1));
			resp = (a_real *)fpu.foo;
			storeres = 1;
			break;
	    case 3:	__asm__("fldt %0" :  : "m" (fpu.st[0])); /* op */
	    		__asm__("fstpl %0" : "=m" (op1));
			op2 = (a_real)(*((double *)fpu.foo));
			loadres = 1;
	    		break;
	}

	/* compute default result (exceptions are masked now) */
	switch (type)
	{   case 1:
	    case 2:	result = op1;
			break;
	    case 3:	if (status & EXC_NXA)	/* inexact - result in st(0) */
			    result = op1;
			else
			    switch (code & 0x0738)
			    {   case 0x400: result = op1+op2; break;
				case 0x420: result = op1-op2; break;
				case 0x408: result = op1*op2; break;
				case 0x430: result = op1/op2; break;
			    }
			break;
	}

        /* test status flag for exception occurrence */
	/* set pxsc flags or generate pxsc exception */
	if (status & EXC_NVA) {
	    if (e_io_e())
		e_trap(INV_OP+E_IEEE, 6, E_TMSG, 69,
					 E_TDBL+E_TEXT(3), &result,
					 E_TDBL|E_TRES, &result);
	    else
		e_sioo();
	}
	
	if (status & EXC_OFA) {
	    if (e_of_e())
		e_trap(OVERFLOW+E_IEEE, 4, E_TDBL+E_TEXT(3), &result,
					   E_TDBL|E_TRES, &result);
	    else
		e_sofo();
	}
	
	if (status & EXC_UFA) {
	    if (e_uf_e())
		e_trap(UNDERFLOW+E_IEEE, 4, E_TDBL+E_TEXT(3), &result,
					    E_TDBL|E_TRES, &result);
	    else
		e_sufo();
	}
	
	if (status & EXC_DZA) {
	    if (e_dz_e())
		e_trap(DIV_BY_ZERO+E_IEEE, 4, E_TDBL+E_TEXT(3), &result,
					      E_TDBL|E_TRES, &result);
	    else
		e_sdzo();
	}
	
	if (status & EXC_NXA) {
	    if (e_ie_e())
		e_trap(INEXACT, 0);
	    else
		e_sieo();
	}

	if (storeres) *resp = result;
        __asm__("frstor %0" :  : "m" (fpu));
	if (loadres) {
	    __asm__("fldcw %0" : "=m" (control));   /* all exceptions masked */
	    __asm__("fldl %0": :  "m" (result));    /* load result */
	    control = fpu.cwd & 0x0f3f;		    /* restore exceptions */
	    __asm__("fldcw %0" : "=m" (control));
	}
        signal(SIGFPE, r_ftrp);		/* restore signal handler */
}

#endif /* IBM_LINUX_C+IBM_EMX_C */

#if IBM_RS6000_C

#define EXC_NVA         0x20000000      /* invalid operation */
#define EXC_OFA         0x10000000      /* overflow */
#define EXC_UFA         0x08000000      /* underflow */
#define EXC_DZA         0x04000000      /* underflow */
#define EXC_NXA         0x02000000      /* inexact */

#ifdef LINT_ARGS
void r_ftrp(int sig, int dummy, struct sigcontext *scp)
#else
void r_ftrp(sig, dummy, scp) 
int sig;
int dummy;
struct sigcontext *scp ;
#endif
{
	int code;			/* exception code */
	struct mstsave *state;		/* context of the interrupted routine */
	fp_sh_info_t flt_context;	/* FPU context, from fp_sh_info() */
	unsigned long op;		/* FPU opcode of trapped instruction */
	double op1, op2, result;	/* operands and result */
	int msg;			/* Meldung bei invalid operation */
	double ieee2;

	state = &scp->sc_jmpbuf.jmp_context;
	fp_sh_info(scp, &flt_context, FP_SH_INFO_SIZE);
	code = flt_context.trap;

	fprintf(stderr, "ieee exception code %x occured at pc %x\n",
		code, state->iar);

	/* determine opcode of failed instruction */
	/* fc21102a add, fc211028 sub, fc2100b2 mul, fc211024 div */
	op = *((unsigned long *)state->iar);

	/* compute correction factor 2**1536 = ieee2*ieee2 */
	*(((unsigned int *)&ieee2)  ) = (768+1023) << 20;
	*(((unsigned int *)&ieee2)+1) = 0;

	/* get operands and determine default result */
	op1 = state->fpr[1];
	op2 = state->fpr[2];
	fp_disable_all();		/* disable all traps */
	switch (flt_context.fpscr & 3)	/* set rounding mode */
	{
	    case 0 : (void)fp_swap_rnd(FP_RND_RN); break;
	    case 1 : (void)fp_swap_rnd(FP_RND_RZ); break;
	    case 2 : (void)fp_swap_rnd(FP_RND_RP); break;
	    case 3 : (void)fp_swap_rnd(FP_RND_RM); break;
	}
	if (!(code & (EXC_NVA|EXC_OFA|EXC_UFA|EXC_DZA)))
	    result = op1;	/* inexact - operation has been performed */
	else if (code & EXC_OFA)	/* overflow - correct result */
	    result = (op1*ieee2)*ieee2;
	else if (code & EXC_UFA)	/* underflow - correct result */
	    result = (op1/ieee2)/ieee2;
	else
	    switch (op)
	    {
		case 0xfc21102a : result = op1+op2; break;
		case 0xfc211028 : result = op1-op2; break;
		case 0xfc2100b2 : result = op1*op2; break;
		case 0xfc211024 : result = op1/op2; break;
		default : fprintf(stderr, "unknown floating-point op code\n");
			    e_trap(INV_OP+E_EXIT, 0);
			    break;
	    }

	/* test status flags for exception occurrence */
	/* set pxsc flags or generate pxsc exception  */
	if (code & EXC_NVA) {
	    msg = 69;				/* keine Meldung */
	    if (code & 0x01000000) msg = 5;	/* SNaN as operand */
	    if (code & 0x00800000) msg = 9;	/* inf - inf */
	    if (code & 0x00400000) msg = 4;	/* inf / inf */
	    if (code & 0x00200000) msg = 2;	/* 0 / 0 */
	    if (code & 0x00100000) msg = 10;	/* inf * 0 */
	    if (e_io_e())
		e_trap(INV_OP+E_IEEE, 8, E_TMSG, msg,
			E_TDBL+E_TEXT(1), &op1, E_TDBL+E_TEXT(2), &op2,
			E_TDBL+E_TEXT(3), &result);
	    else
		e_sioo();
	}
	
	if (code & EXC_OFA) {
	    if (e_of_e())
		e_trap(OVERFLOW+E_IEEE, 6,
			E_TDBL+E_TEXT(1), &op1, E_TDBL+E_TEXT(2), &op2,
			E_TDBL+E_TEXT(3), &result);
	    else
		e_sofo();
	}
	
	if (code & EXC_UFA) {
	    if (e_uf_e())
		e_trap(UNDERFLOW+E_IEEE, 6,
			E_TDBL+E_TEXT(1), &op1, E_TDBL+E_TEXT(2), &op2,
			E_TDBL+E_TEXT(3), &result);
	    else
		e_sufo();
	}
	
	if (code & EXC_DZA) {
	    if (e_dz_e())
		e_trap(DIV_BY_ZERO+E_IEEE, 6,
			E_TDBL+E_TEXT(1), &op1, E_TDBL+E_TEXT(2), &op2,
			E_TDBL+E_TEXT(3), &result);
	    else
		e_sdzo();
	}
	
	if (code & EXC_NXA) {
	    if (e_ie_e())
		e_trap(INEXACT, 0);
	    else
		e_sieo();
	}

	/* store result */
	state->fpr[1] = result;

	/* clear exception bits to prevent recurrence of the trap */
	fp_sh_set_stat(scp,
	    (flt_context.fpscr & ((fpstat_t) ~flt_context.trap)));

	/* increment return address if it points to trapped instruction */
	if (flt_context.flags & FP_IAR_STAT)
	    state->iar += 4;
	return;
}

#endif /* IBM_RS6000_C */

#else     /* IEEE_HARDWARE */
/* dummy entry */
#ifdef LINT_ARGS
void r_ftrp(void)
#else
void r_ftrp()
#endif
{
}
#endif	/* IEEE_HARDWARE */





