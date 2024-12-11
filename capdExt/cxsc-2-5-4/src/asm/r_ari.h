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

/* CVS $Id: r_ari.h,v 1.12 2014/01/30 17:23:52 cxsc Exp $ */

#ifndef _CXSC_R_ARI_H_INCLUDED
#define _CXSC_R_ARI_H_INCLUDED

namespace cxsc {


#if defined(CXSC_x86_64) || defined(CXSC_x86_i)

static u_int32_t cxscup = 0x5F80;
static u_int32_t cxscdown = 0x3F80;
static u_int32_t cxscnear = 0x1F80;

inline double ra_addu(double a,double b) {
   #ifdef CXSC_RND_STORE
     u_int32_t store = cxscnear; 
     asm volatile ("stmxcsr %4\n" "ldmxcsr %2\n" "addsd %1, %0\n" "ldmxcsr %3\n" : "+x" (a) : "x" (b), "m" (cxscup), "m" (store), "m" (store) );
   #else
     asm volatile ("ldmxcsr %2\n" "addsd %1, %0\n" "ldmxcsr %3\n" : "+x" (a) : "x" (b), "m" (cxscup), "m" (cxscnear) );
   #endif
   return a;
}

inline double ra_addd(double a,double b) {
   #ifdef CXSC_RND_STORE
     u_int32_t store = cxscnear; 
     asm volatile ("stmxcsr %4\n" "ldmxcsr %2\n" "addsd %1, %0\n" "ldmxcsr %3\n" : "+x" (a) : "x" (b), "m" (cxscdown), "m" (store), "m" (store) );
   #else
     asm volatile ("ldmxcsr %2\n" "addsd %1, %0\n" "ldmxcsr %3\n" : "+x" (a) : "x" (b), "m" (cxscdown), "m" (cxscnear) );
   #endif
   return a;
}

inline double ra_subu(double a,double b) {
   #ifdef CXSC_RND_STORE
     u_int32_t store = cxscnear; 
     asm volatile ("stmxcsr %4\n" "ldmxcsr %2\n" "subsd %1, %0\n" "ldmxcsr %3\n" : "+x" (a) : "x" (b), "m" (cxscup), "m" (store), "m" (store) );
   #else
     asm volatile ("ldmxcsr %2\n" "subsd %1, %0\n" "ldmxcsr %3\n" : "+x" (a) : "x" (b), "m" (cxscup), "m" (cxscnear) );
   #endif
   return a;
}

inline double ra_subd(double a,double b) {
   #ifdef CXSC_RND_STORE
     u_int32_t store = cxscnear; 
     asm volatile ("stmxcsr %4\n" "ldmxcsr %2\n" "subsd %1, %0\n" "ldmxcsr %3\n" : "+x" (a) : "x" (b), "m" (cxscdown), "m" (store), "m" (store) );
   #else
     asm volatile ("ldmxcsr %2\n" "subsd %1, %0\n" "ldmxcsr %3\n" : "+x" (a) : "x" (b), "m" (cxscdown), "m" (cxscnear) );
   #endif
   return a;
}

inline double ra_mulu(double a,double b) {
   #ifdef CXSC_RND_STORE
     u_int32_t store = cxscnear; 
     asm volatile ("stmxcsr %4\n" "ldmxcsr %2\n" "mulsd %1, %0\n" "ldmxcsr %3\n" : "+x" (a) : "x" (b), "m" (cxscup), "m" (store), "m" (store) );
   #else
     asm volatile ("ldmxcsr %2\n" "mulsd %1, %0\n" "ldmxcsr %3\n" : "+x" (a) : "x" (b), "m" (cxscup), "m" (cxscnear) );
   #endif
   return a;
}

inline double ra_muld(double a,double b) {
   #ifdef CXSC_RND_STORE
     u_int32_t store = cxscnear; 
     asm volatile ("stmxcsr %4\n" "ldmxcsr %2\n" "mulsd %1, %0\n" "ldmxcsr %3\n" : "+x" (a) : "x" (b), "m" (cxscdown), "m" (store), "m" (store) );
   #else
     asm volatile ("ldmxcsr %2\n" "mulsd %1, %0\n" "ldmxcsr %3\n" : "+x" (a) : "x" (b), "m" (cxscdown), "m" (cxscnear) );
   #endif
   return a;
}

inline double ra_divu(double a,double b) {
   #ifdef CXSC_RND_STORE
     u_int32_t store = cxscnear; 
     asm volatile ("stmxcsr %4\n" "ldmxcsr %2\n" "divsd %1, %0\n" "ldmxcsr %3\n" : "+x" (a) : "x" (b), "m" (cxscup), "m" (store), "m" (store) );
   #else
     asm volatile ("ldmxcsr %2\n" "divsd %1, %0\n" "ldmxcsr %3\n" : "+x" (a) : "x" (b), "m" (cxscup), "m" (cxscnear) );
   #endif
   return a;
}

inline double ra_divd(double a,double b) {
   #ifdef CXSC_RND_STORE
     u_int32_t store = cxscnear; 
     asm volatile ("stmxcsr %4\n" "ldmxcsr %2\n" "divsd %1, %0\n" "ldmxcsr %3\n" : "+x" (a) : "x" (b), "m" (cxscdown), "m" (store), "m" (store) );
   #else
     asm volatile ("ldmxcsr %2\n" "divsd %1, %0\n" "ldmxcsr %3\n" : "+x" (a) : "x" (b), "m" (cxscdown), "m" (cxscnear) );
   #endif
   return a;
}

#elif defined(CXSC_x86)

static u_int16_t cxscup = 0xB3F;
static u_int16_t cxscdown = 0x73F;
static u_int16_t cxscnear = 0x33F;

inline double ra_addu(double a,double b) {
   double c;
   #ifdef CXSC_RND_STORE
     u_int16_t store = cxscnear; 
     asm volatile ("fstcw %4\n" "fldcw %3\n" "fadd %2, %1\n" "fstpl %0\n" "fldcw %5\n" : "=m" (c) : "t" (a), "u" (b), "m" (cxscup), "m" (store), "m" (store) : "st");
   #else
     asm volatile ("fldcw %3\n" "fadd %2, %1\n" "fstpl %0\n" "fldcw %4\n" : "=m" (c) : "t" (a), "u" (b), "m" (cxscup), "m" (cxscnear)  : "st");
   #endif
   return c;
}

inline double ra_addd(double a,double b) {
   double c;
   #ifdef CXSC_RND_STORE
     u_int16_t store = cxscnear; 
     asm volatile ("fstcw %4\n" "fldcw %3\n" "fadd %2, %1\n" "fstpl %0\n" "fldcw %5\n" : "=m" (c) : "t" (a), "u" (b), "m" (cxscdown), "m" (store), "m" (store) : "st");
   #else
     asm volatile ("fldcw %3\n" "fadd %2, %1\n" "fstpl %0\n" "fldcw %4\n" : "=m" (c) : "t" (a), "u" (b), "m" (cxscdown), "m" (cxscnear)  : "st");
   #endif
   return c;
}

inline double ra_subu(double a,double b) {
   double c;
   #ifdef CXSC_RND_STORE
     u_int16_t store = cxscnear; 
     asm volatile ("fstcw %4\n" "fldcw %3\n" "fsub %2, %1\n" "fstpl %0\n" "fldcw %5\n" : "=m" (c) : "t" (a), "u" (b), "m" (cxscup), "m" (store), "m" (store) : "st");
   #else
     asm volatile ("fldcw %3\n" "fsub %2, %1\n" "fstpl %0\n" "fldcw %4\n" : "=m" (c) : "t" (a), "u" (b), "m" (cxscup), "m" (cxscnear)  : "st");
   #endif
   return c;
}

inline double ra_subd(double a,double b) {
   double c;
   #ifdef CXSC_RND_STORE
     u_int16_t store = cxscnear; 
     asm volatile ("fstcw %4\n" "fldcw %3\n" "fsub %2, %1\n" "fstpl %0\n" "fldcw %5\n" : "=m" (c) : "t" (a), "u" (b), "m" (cxscdown), "m" (store), "m" (store) : "st");
   #else
     asm volatile ("fldcw %3\n" "fsub %2, %1\n" "fstpl %0\n" "fldcw %4\n" : "=m" (c) : "t" (a), "u" (b), "m" (cxscdown), "m" (cxscnear)  : "st");
   #endif
   return c;
}

inline double ra_mulu(double a,double b) {
   double c;
   #ifdef CXSC_RND_STORE
     u_int16_t store = cxscnear; 
     asm volatile ("fstcw %4\n" "fldcw %3\n" "fmul %2, %1\n" "fstpl %0\n" "fldcw %5\n" : "=m" (c) : "t" (a), "u" (b), "m" (cxscup), "m" (store), "m" (store) : "st");
   #else
     asm volatile ("fldcw %3\n" "fmul %2, %1\n" "fstpl %0\n" "fldcw %4\n" : "=m" (c) : "t" (a), "u" (b), "m" (cxscup), "m" (cxscnear)  : "st");
   #endif
   return c;
}

inline double ra_muld(double a,double b) {
   double c;
   #ifdef CXSC_RND_STORE
     u_int16_t store = cxscnear; 
     asm volatile ("fstcw %4\n" "fldcw %3\n" "fmul %2, %1\n" "fstpl %0\n" "fldcw %5\n" : "=m" (c) : "t" (a), "u" (b), "m" (cxscdown), "m" (store), "m" (store) : "st");
   #else
     asm volatile ("fldcw %3\n" "fmul %2, %1\n" "fstpl %0\n" "fldcw %4\n" : "=m" (c) : "t" (a), "u" (b), "m" (cxscdown), "m" (cxscnear)  : "st");
   #endif
   return c;
}

inline double ra_divu(double a,double b) {
   double c;
   #ifdef CXSC_RND_STORE
     u_int16_t store = cxscnear; 
     asm volatile ("fstcw %4\n" "fldcw %3\n" "fdiv %2, %1\n" "fstpl %0\n" "fldcw %5\n" : "=m" (c) : "t" (a), "u" (b), "m" (cxscup), "m" (store), "m" (store) : "st");
   #else
     asm volatile ("fldcw %3\n" "fdiv %2, %1\n" "fstpl %0\n" "fldcw %4\n" : "=m" (c) : "t" (a), "u" (b), "m" (cxscup), "m" (cxscnear)  : "st");
   #endif
   return c;
}

inline double ra_divd(double a,double b) {
   double c;
   #ifdef CXSC_RND_STORE
     u_int16_t store = cxscnear; 
     asm volatile ("fstcw %4\n" "fldcw %3\n" "fdiv %2, %1\n" "fstpl %0\n" "fldcw %5\n" : "=m" (c) : "t" (a), "u" (b), "m" (cxscdown), "m" (store), "m" (store) : "st");
   #else
     asm volatile ("fldcw %3\n" "fdiv %2, %1\n" "fstpl %0\n" "fldcw %4\n" : "=m" (c) : "t" (a), "u" (b), "m" (cxscdown), "m" (cxscnear)  : "st");
   #endif
   return c;
}

#elif defined(CXSC_PPC64)

inline double ra_addu(double a,double b) {
   double r(0);
   #ifdef CXSC_RND_STORE
     double store;
     asm volatile ("mffs %3\n" "mtfsfi 7,2\n" "fadd %0, %1, %2\n" "mtfsf 0x03,%3\n" : "=f" (r) : "f"(a), "f" (b), "f"(store) );
   #else
     asm volatile ("mtfsfi 7,2\n" "fadd %0, %1, %2\n" "mtfsfi 7,0\n" : "=f" (r) : "f"(a), "f" (b) );
   #endif
   return r;
}

inline double ra_addd(double a,double b) {
   double r(0);
   #ifdef CXSC_RND_STORE
     double store;
     asm volatile ("mffs %3\n" "mtfsfi 7,3\n" "fadd %0, %1, %2\n" "mtfsf 0x03,%3\n" : "=f" (r) : "f"(a), "f" (b), "f"(store) );
   #else
     asm volatile ("mtfsfi 7,3\n" "fadd %0, %1, %2\n" "mtfsfi 7,0\n" : "=f" (r) : "f"(a), "f" (b) );
   #endif
   return r;
}

inline double ra_subu(double a,double b) {
   double r(0);
   #ifdef CXSC_RND_STORE
     double store;
     asm volatile ("mffs %3\n" "mtfsfi 7,2\n" "fsub %0, %1, %2\n" "mtfsf 0x03,%3\n" : "=f" (r) : "f"(a), "f" (b), "f"(store) );
   #else
     asm volatile ("mtfsfi 7,2\n" "fsub %0, %1, %2\n" "mtfsfi 7,0\n" : "=f" (r) : "f"(a), "f" (b) );
   #endif
   return r;
}

inline double ra_subd(double a,double b) {
   double r(0);
   #ifdef CXSC_RND_STORE
     double store;
     asm volatile ("mffs %3\n" "mtfsfi 7,3\n" "fsub %0, %1, %2\n" "mtfsf 0x03,%3\n" : "=f" (r) : "f"(a), "f" (b), "f"(store) );
   #else
     asm volatile ("mtfsfi 7,3\n" "fsub %0, %1, %2\n" "mtfsfi 7,0\n" : "=f" (r) : "f"(a), "f" (b) );
   #endif
   return r;
}

inline double ra_mulu(double a,double b) {
   double r(0);
   #ifdef CXSC_RND_STORE
     double store;
     asm volatile ("mffs %3\n" "mtfsfi 7,2\n" "fmul %0, %1, %2\n" "mtfsf 0x03,%3\n" : "=f" (r) : "f"(a), "f" (b), "f"(store) );
   #else
     asm volatile ("mtfsfi 7,2\n" "fmul %0, %1, %2\n" "mtfsfi 7,0\n" : "=f" (r) : "f"(a), "f" (b) );
   #endif
   return r;
}

inline double ra_muld(double a,double b) {
   double r(0);
   #ifdef CXSC_RND_STORE
     double store;
     asm volatile ("mffs %3\n" "mtfsfi 7,3\n" "fmul %0, %1, %2\n" "mtfsf 0x03,%3\n" : "=f" (r) : "f"(a), "f" (b), "f"(store) );
   #else
     asm volatile ("mtfsfi 7,3\n" "fmul %0, %1, %2\n" "mtfsfi 7,0\n" : "=f" (r) : "f"(a), "f" (b) );
   #endif
   return r;
}

inline double ra_divu(double a,double b) {
   double r(0);
   #ifdef CXSC_RND_STORE
     double store;
     asm volatile ("mffs %3\n" "mtfsfi 7,2\n" "fdiv %0, %1, %2\n" "mtfsf 0x03,%3\n" : "=f" (r) : "f"(a), "f" (b), "f"(store) );
   #else
     asm volatile ("mtfsfi 7,2\n" "fdiv %0, %1, %2\n" "mtfsfi 7,0\n" : "=f" (r) : "f"(a), "f" (b) );
   #endif
   return r;
}

inline double ra_divd(double a,double b) {
   double r(0);
   #ifdef CXSC_RND_STORE
     double store;
     asm volatile ("mffs %3\n" "mtfsfi 7,3\n" "fdiv %0, %1, %2\n" "mtfsf 0x03,%3\n" : "=f" (r) : "f"(a), "f" (b), "f"(store) );
   #else
     asm volatile ("mtfsfi 7,3\n" "fdiv %0, %1, %2\n" "mtfsfi 7,0\n" : "=f" (r) : "f"(a), "f" (b) );
   #endif
   return r;
}


#endif

}


#endif

