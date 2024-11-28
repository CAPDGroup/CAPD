# -*- Makefile -*-

AM_CPPFLAGS = -I${top_srcdir}/src/asm

noinst_LTLIBRARIES = librts.la

librts_la_SOURCES = \
	a_abs_.c a_add_.c a_band.c a_bclr.c a_beqv.c a_bmsb.c a_bnot.c \
	a_bons.c a_bor_.c a_brtt.c a_bset.c a_bshf.c a_btst.c a_bxor.c \
	a_clck.c a_div_.c a_exit.c a_free.c a_gets.c a_ival.c a_ixch.c \
	a_lloc.c a_mod_.c a_mul_.c a_nilc.c a_popt.c a_pptr.c a_sets.c \
	a_sqr_.c a_sub_.c a_sval.c a_syst.c a_tick.c a_time.c a_umin.c \
	b_accu.c b_acsh.c b_acth.c b_addc.c b_addm.c b_addu.c b_adj.c \
	b_adpp.c b_asgn.c b_asin.c b_asiv.c b_asnh.c b_atan.c b_atav.c \
	b_atnh.c b_baad.c b_bacm.c b_badd.c b_badj.c b_ball.c b_banx.c \
	b_basu.c b_bcad.c b_bcat.c b_bcdi.c b_bcid.c b_bclr.c b_bcmp.c \
	b_bcpy.c b_bcsu.c b_bdiv.c b_bdvn.c b_bini.c b_biv2.c b_biv_.c \
	b_bivp.c b_bldx.c b_blgx.c b_bmat.c b_bmcm.c b_bmdv.c b_bms1.c \
	b_bms2.c b_bmsh.c b_bmsp.c b_bmts.c b_bmul.c b_bmun.c b_bnxt.c \
	b_bpnt.c b_brnd.c b_bshf.c b_bsub.c b_busp.c b_chck.c b_cmcp.c \
	b_comp.c b_conf.c b_coni.c b_cos_.c b_cosh.c b_cot_.c b_coth.c \
	b_cscn.c b_deko.c b_dtol.c b_errr.c b_exp_.c b_expe.c b_form.c \
	b_geta.c b_gini.c b_glbl.c b_heap.c b_ifrm.c b_inv1.c b_inv2.c \
	b_irnd.c b_kapp.c b_lnva.c b_lnve.c b_log_.c b_loga.c b_ltod.c \
	b_ltor.c b_mdiv.c b_muad.c b_op88.c b_out.c b_outf.c b_outi.c \
	b_outm.c b_pi__.c b_pign.c b_pool.c b_popt.c b_pow_.c b_prod.c \
	b_rnd.c b_rndd.c b_rndn.c b_rndu.c b_rtol.c b_scan.c b_shl1.c \
	b_shlu.c b_shr1.c b_shru.c b_sico.c b_sin_.c b_sinh.c b_snhv.c \
	b_sqrt.c b_sqrv.c b_subc.c b_subm.c b_subu.c b_tadd.c b_tadj.c \
	b_tan_.c b_tanh.c b_tcom.c b_tdek.c b_tdiv.c b_test.c b_text.c \
	b_tmdv.c b_tmpf.c b_tmul.c b_trnd.c b_tsub.c c_cadd.c c_csub.c \
	c_dadd.c c_dsta.c c_dsub.c c_padd.c c_psub.c c_rcad.c c_rcsb.c \
	c_scps.c c_scpy.c c_stad.c c_stan.c c_stau.c d_ass.c d_assc.c \
	d_assi.c d_assz.c d_chs.c d_clr.c d_dadd.c d_dsub.c d_eq.c d_free.c \
	d_ge.c d_gt.c d_init.c d_le.c d_lt.c d_ne.c d_padd.c d_psub.c \
	d_radd.c d_rsub.c d_sign.c d_stad.c d_stan.c d_stau.c d_temp.c \
	d_utmp.c d_vlcp.c e_actn.c e_back.c e_bmsg.c e_data.c e_ieee.c \
	e_ienv.c e_popp.c e_push.c e_tmrt.c e_tmsg.c e_tprt.c e_trap.c \
	e_xall.c e_xarg.c e_xdbz.c e_xine.c e_xiob.c e_xiop.c e_xios.c \
	e_xnor.c e_xofl.c e_xset.c e_xufl.c f_args.c f_assg.c f_back.c \
	f_bhex.c f_eofp.c f_exst.c f_flsh.c f_free.c f_get_.c f_getc.c \
	f_op88.c f_popt.c f_put_.c f_putc.c f_quer.c f_rdf1.c f_rdh1.c \
	f_rdi1.c f_rdl1.c f_rdl2.c f_rdln.c f_rdr1.c f_rdr2.c f_rds1.c \
	f_rdv1.c f_read.c f_rhex.c f_rint.c f_rset.c f_rstn.c f_rwri.c \
	f_rwrn.c f_sexs.c f_srse.c f_srwi.c f_whex.c f_wint.c f_wrb1.c \
	f_wrb2.c f_wrc1.c f_wrc2.c f_wrd1.c f_wrf1.c f_wrf2.c f_wrf3.c \
	f_wrf4.c f_wrh1.c f_wrh2.c f_wri1.c f_wri2.c f_wrid.c f_writ.c \
	f_wrl1.c f_wrl2.c f_wrl3.c f_wrl4.c f_wrln.c f_wrr1.c f_wrr2.c \
	f_wrr3.c f_wrr4.c f_wrs1.c f_wrs2.c i_acos.c i_acot.c i_acsh.c \
	i_acth.c i_asin.c i_asnh.c i_atan.c i_atnh.c i_cns2.c i_cnst.c \
	i_cos.c i_cosh.c i_cot.c i_coth.c i_dadd.c i_dsta.c i_dsub.c i_ep10.c \
	i_exp.c i_exp2.c i_iadd.c i_ista.c i_isub.c i_lg10.c i_log.c i_log2.c \
	i_loga.c i_padd.c i_pow.c i_psub.c i_read.c i_riad.c i_risb.c \
	i_scps.c i_scpy.c i_sin.c i_sinh.c i_sqrt.c i_tan.c i_tanh.c l_abs.c \
	l_acos.c l_acot.c l_acsh.c l_acth.c l_addc.c l_addd.c l_addu.c \
	l_asin.c l_asnh.c l_ass.c l_atan.c l_atn2.c l_atnh.c l_comp.c \
	l_conv.c l_cos.c l_cosh.c l_cot.c l_coth.c l_divc.c l_divd.c l_divu.c \
	l_eq.c l_exct.c l_exp.c l_expo.c l_flot.c l_free.c l_ge.c l_gt.c \
	l_init.c l_le.c l_log.c l_loga.c l_lt.c l_mant.c l_mlen.c l_mulc.c \
	l_muld.c l_mulu.c l_ne.c l_pow.c l_prec.c l_pred.c l_read.c l_rond.c \
	l_rprc.c l_rval.c l_sign.c l_sin.c l_sinh.c l_sqrt.c l_subc.c \
	l_subd.c l_subu.c l_succ.c l_tan.c l_tanh.c l_temp.c l_trun.c \
	l_ttmp.c l_umin.c l_utmp.c l_vlcp.c l_whex.c l_writ.c o_user.c \
	p_init.c r_abs.c r_acos.c r_acot.c r_acsh.c r_acth.c r_addd.c \
	r_addn.c r_addu.c r_asin.c r_asnh.c r_ass.c r_atan.c r_atn2.c \
	r_atnh.c r_aval.c r_ceil.c r_clss.c r_cnsd.c r_cnst.c r_cnsu.c \
	r_comp.c r_conv.c r_cos.c r_cosh.c r_cot.c r_coth.c r_divd.c r_divn.c \
	r_divu.c r_ep10.c r_eq.c r_exp.c r_exp2.c r_expo.c r_fini.c r_flor.c \
	r_flot.c r_frac.c r_ftrp.c r_ge.c r_gt.c r_le.c r_lg10.c r_log.c \
	r_log2.c r_loga.c r_lt.c r_mant.c r_muld.c r_muln.c r_mulu.c r_ne.c \
	r_outp.c r_pcmp.c r_pow.c r_pred.c r_rdcr.c r_read.c r_rond.c \
	r_rval.c r_scps.c r_scpy.c r_sign.c r_sin.c r_sinh.c r_sqr_.c \
	r_sqrt.c r_subd.c r_subn.c r_subu.c r_succ.c r_sval.c r_tan.c \
	r_tanh.c r_trun.c r_umin.c r_valu.c r_writ.c s_aacc.c s_aaeq.c \
	s_aage.c s_aagt.c s_aain.c s_aale.c s_aalt.c s_aane.c s_accc.c \
	s_aceq.c s_acge.c s_acgt.c s_acle.c s_aclt.c s_acne.c s_add.c \
	s_aimg.c s_alen.c s_asgn.c s_cacc.c s_caeq.c s_cage.c s_cagt.c \
	s_cain.c s_cale.c s_calt.c s_cane.c s_char.c s_conc.c s_cons.c \
	s_cpos.c s_cseq.c s_csge.c s_csgt.c s_csle.c s_cslt.c s_csne.c \
	s_date.c s_eteq.c s_etge.c s_etgt.c s_etin.c s_etle.c s_etlt.c \
	s_etne.c s_free.c s_genv.c s_init.c s_ins1.c s_ins2.c s_int_.c \
	s_inxc.c s_ixch.c s_ixcn.c s_len_.c s_mul.c s_real.c s_rhex.c \
	s_sceq.c s_scge.c s_scgt.c s_scle.c s_sclt.c s_scne.c s_slen.c \
	s_spos.c s_sseq.c s_ssge.c s_ssgt.c s_ssle.c s_sslt.c s_ssne.c \
	s_stat.c s_sub.c s_suba.c s_subs.c s_temp.c s_utmp.c s_vlcp.c \
	s_whex.c s_zero.c t_abse.c t_acos.c t_acot.c t_acsh.c t_acth.c \
	t_acvt.c t_asin.c t_asnh.c t_asvt.c t_atan.c t_atn2.c t_atnh.c \
	t_cha1.c t_cha2.c t_cher.c t_chi1.c t_chi2.c t_chiv.c t_chk1.c \
	t_chse.c t_cmpe.c t_cnst.c t_conv.c t_cos.c t_cosh.c t_cot.c t_coth.c \
	t_cpye.c t_ctoe.c t_dadd.c t_das_.c t_dcmp.c t_dcnv.c t_dcpy.c \
	t_dcst.c t_deml.c t_dhi1.c t_dhi2.c t_dhl1.c t_dhl2.c t_dmsc.c \
	t_dmul.c t_dnrm.c t_dsub.c t_ecst.c t_ehi1.c t_ehi2.c t_ehl1.c \
	t_ehl2.c t_etoa.c t_etoi.c t_ex10.c t_exme.c t_exp.c t_exp2.c \
	t_glbl.c t_gtrg.c t_iaah.c t_iabs.c t_iabt.c t_iacs.c t_iact.c \
	t_ianh.c t_iash.c t_iasn.c t_iath.c t_iatn.c t_ichs.c t_icos.c \
	t_icot.c t_icpe.c t_icsh.c t_icth.c t_ie10.c t_iex2.c t_iexp.c \
	t_il10.c t_ilg2.c t_ilog.c t_ipow.c t_irnd.c t_is_c.c t_isin.c \
	t_isnh.c t_isqt.c t_itan.c t_itnh.c t_lg10.c t_lnpe.c t_log.c \
	t_log2.c t_md2e.c t_merr.c t_mtyp.c t_mupi.c t_pow.c t_rtrg.c \
	t_s2xm.c t_s_ei.c t_s_ln.c t_s_lp.c t_satn.c t_scee.c t_sico.c \
	t_sin.c t_sinh.c t_sinv.c t_sle_.c t_sqme.c t_sqrt.c t_sqt_.c \
	t_srnd.c t_ssin.c t_stan.c t_tan.c t_tanh.c t_tnct.c t_xtpe.c \
	t_xtre.c y_alck.c y_asgn.c y_free.c y_inid.c y_init.c y_inxn.c \
	y_ixch.c y_ixcn.c y_new.c y_stat.c y_suba.c y_temp.c y_utmp.c \
	y_vlcp.c y_ynxn.c y_yxch.c y_yxcn.c z_ciad.c z_cisb.c z_czad.c \
	z_czsb.c z_dadd.c z_dsta.c z_dsub.c z_izad.c z_izsb.c z_padd.c \
	z_psub.c z_rzad.c z_rzsb.c z_scps.c z_scpy.c z_zadd.c z_zsta.c \
	z_zsub.c
