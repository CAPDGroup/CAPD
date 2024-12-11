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

/* CVS $Id: o_name.h,v 1.21 2014/01/30 17:24:10 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : o_name.h                              */
/*                                                              */
/*                   fsrqn,fsrqu,fsrqd,fsrqc,temd =b=           */
/****************************************************************/

#define addm            b_addm
#define adduc           b_addu
#define adjust          b_adj
#define addii           b_badd
#define allocu          b_ball
#define anext           b_banx
#define asub            b_basu
#define cadd            b_bcad
#define caddt           b_bcat
#define convdi          b_bcdi
#define convid          b_bcid
#define cleari          b_bclr
#define compii          b_bcmp
#define copyii          b_bcpy
#define csub            b_bcsu
#define divii           b_bdiv
#define divin           b_bdvn
#define initi           b_bini
#define i_iv            b_biv_
#define i_inv2          b_biv2
#define i_invp          b_bivp
#define i_ldex          b_bldx
#define i_lgex          b_blgx
#define madd            b_bmad
#define mapiu           b_bmap
#define maddt           b_bmat
#define madd11          b_bma1
#define madd12          b_bma2
#define mcomp           b_bmcm
#define mcopy           b_bmcp
#define mdiv32          b_bmdv
#define mshift          b_bmsh
#define mspiu           b_bmsp
#define msub            b_bmsu
#define msub11          b_bms1
#define msub12          b_bms2
#define mtest           b_bmts
#define mulii           b_bmul
#define mulin           b_bmun
#define nexti           b_bnxt
#define i_point         b_bpnt
#define roundi          b_brnd
#define shifti          b_bshf
#define subii           b_bsub
#define mapuu           b_buap
#define clearu          b_bucl
#define mspuu           b_busp
/*
#define check           b_chck
*/
#define compos          b_comp
#define conf            b_conf
#define coni            b_coni
#define deko            b_deko
#define r_form          b_form
#define getaccu         b_geta
#define invoke          b_inv1
#define invoke2         b_inv2
#define l_to_r          b_ltor
#define mantdiv         b_mdiv
#define muladd          b_muad
#define out             b_out
#define outi            b_outi
#define round           b_rnd
#define roundd          b_rndd
#define roundn          b_rndn
#define roundu          b_rndu
#define r_to_l          b_rtol
#define rfscan          b_scan
#define r_scan          b_sscn
#define shftul          b_shlu
#define shiftl          b_shl1
#define shftur          b_shru
#define shiftr          b_shr1
#define subc            b_subc
#define subm            b_subm
#define subuc           b_subu
#define testu           b_test

#define L4ovPi          b_4opi
#define Larccos         b_acos
#define Larccot         b_acot
#define Larcosh         b_acsh
#define Larcoth         b_acth
#define Lassign         b_asgn
#define Larcsin         b_asin
#define arcsinve        b_asiv
#define Larsinh         b_asnh
#define Larctan         b_atan
#define arctanve        b_atav
#define Latan2          b_atn2
#define Lartanh         b_atnh
#define Lcos            b_cos_
#define Lcosh           b_cosh
#define Lcot            b_cot_
#define Lcoth           b_coth
#define Lcurrprec       b_cprc
#define Lvardrop        b_drop
#define Lerrmsg         b_errm
#define Lerror          b_errr
#define EUFfac          b_euff
#define Lexp            b_exp_
#define expve           b_expe
#define Faddeps         r_succ
#define Fln2            b_fln2
#define Flnb            b_flnb
#define Fsubeps         r_pred
#define FUFfac          b_fuff
#define FIUFfac         b_fiuf
#define Lvarget         b_get_
#define Lgiflag         b_gifl
#define Lginit          b_gini
#define Ldebug          b_ldbg
#define Leps            b_leps
#define LFminsq         b_lfmi
#define LFulp           b_lflp
#define LMinreal        b_lmrl
#define LhD             b_lhd_
#define LhE             b_lhe_
#define LhF             b_lhf_
#define LhI             b_lhi_
#define LhV             b_lhv_
#define Lintern         b_lint
#define Lmindbl         b_lmin
#define Llnbase         b_lnbs
#define lnvea           b_lnva
#define lnve            b_lnve
#define Lln             b_log_
#define Lloga           b_loga
#define Lone            b_lone
#define Lrnd            b_lrnd
#define Lversion        b_lver
#define Lpi             b_pi__
#define LPiov4          b_pio4
#define LpiGen          b_pign
#define Lpower          b_pow_
#define Lroutine        b_rout
#define sicovea         b_sico
#define Lsin            b_sin_
#define Lsinh           b_sinh
#define sinhvea         b_snhv
#define Lsqrt           b_sqrt
#define sqrtve          b_sqrv
#define Ltan            b_tan_
#define Ltanh           b_tanh

#define cm              b_cm__
#define cp              b_cp__
#define Maxl            b_maxl

#define io_enabled      e_io_e
#define dz_enabled      e_dz_e
#define of_enabled      e_of_e
#define uf_enabled      e_uf_e
#define ie_enabled      e_ie_e

#define io_occurred     e_io_o
#define dz_occurred     e_dz_o
#define of_occurred     e_of_o
#define uf_occurred     e_uf_o
#define ie_occurred     e_ie_o

#define io_set          e_sioo
#define dz_set          e_sdzo
#define of_set          e_sofo
#define uf_set          e_sufo
#define ie_set          e_sieo

#define io_reset        e_rioo
#define dz_reset        e_rdzo
#define of_reset        e_rofo
#define uf_reset        e_rufo
#define ie_reset        e_rieo

#define all_set         e_sall
#define all_reset       e_rall
#define all_save        e_save
#define all_restore     e_rest

#define u_ass           l_ass
#define u_free          l_free
#define u_init          l_init
#define u_eq            l_eq
#define u_ge            l_ge
#define u_gt            l_gt
#define u_le            l_le
#define u_lt            l_lt
#define u_ne            l_ne





