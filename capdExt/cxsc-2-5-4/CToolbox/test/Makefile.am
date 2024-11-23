# -*- Makefile -*-

AM_CPPFLAGS = -I${top_srcdir}/CToolbox/Modules -I${top_srcdir}/src \
              -I${top_srcdir}/src/rts -I${top_srcdir}/src/fi_lib \
              -I${top_srcdir}/src/asm

check_PROGRAMS = sumup tester

LDADD = ../../src/libcxsc.la

sumup_SOURCES = sumup.cpp
tester_SOURCES = tester.cpp

P = ${top_srcdir}/CToolbox/Programs

TESTS = tmain

.PHONY: tmain

tmain: ${check_PROGRAMS}
	${P}/cpz_ex <${P}/cpz_ex.in1 >cpz_ex.out1
	${P}/cpz_ex <${P}/cpz_ex.in2 >cpz_ex.out2
	${P}/gop1_ex <${P}/gop1_ex.in >gop1_ex.out
	${P}/gop_ex <${P}/gop_ex.in >gop_ex.out
	${P}/hess_ex <${P}/hess_ex.in >hess_ex.out
	${P}/jac_ex <${P}/jac_ex.in >jac_ex.out
	${P}/lop_ex <${P}/lop_ex.in1 >lop_ex.out1
	${P}/lop_ex <${P}/lop_ex.in2 >lop_ex.out2
	${P}/lss_ex <${P}/lss_ex.in1 >lss_ex.out1
	${P}/lss_ex <${P}/lss_ex.in2 >lss_ex.out2
	${P}/nlfz_ex <${P}/nlfz_ex.in >nlfz_ex.out
	${P}/nlss_ex <${P}/nlss_ex.in1 >nlss_ex.out1
	${P}/nlss_ex <${P}/nlss_ex.in2 >nlss_ex.out2
	${P}/rpe_ex <${P}/rpe_ex.in1 >rpe_ex.out1
	${P}/rpe_ex <${P}/rpe_ex.in2 >rpe_ex.out2
	${P}/xev_ex1 <${P}/xev_ex1.in >xev_ex1.out
	${P}/xev_ex2 <${P}/xev_ex2.in >xev_ex2.out
	${P}/xev_ex3 <${P}/xev_ex3.in >xev_ex3.out
	rm -f ergebnis.dat
	touch ergebnis.dat
	./tester cpz_ex.out1 ${srcdir}/cpz_ex.vgl1
	./tester cpz_ex.out2 ${srcdir}/cpz_ex.vgl2
	./tester gop1_ex.out ${srcdir}/gop1_ex.vgl
	./tester gop_ex.out ${srcdir}/gop_ex.vgl
	./tester hess_ex.out ${srcdir}/hess_ex.vgl
	./tester jac_ex.out ${srcdir}/jac_ex.vgl
	./tester lop_ex.out1 ${srcdir}/lop_ex.vgl1
	./tester lop_ex.out2 ${srcdir}/lop_ex.vgl2
	./tester lss_ex.out1 ${srcdir}/lss_ex.vgl1
	./tester lss_ex.out2 ${srcdir}/lss_ex.vgl2
	./tester nlfz_ex.out ${srcdir}/nlfz_ex.vgl
	./tester nlss_ex.out1 ${srcdir}/nlss_ex.vgl1
	./tester nlss_ex.out2 ${srcdir}/nlss_ex.vgl2
	./tester rpe_ex.out1 ${srcdir}/rpe_ex.vgl1
	./tester rpe_ex.out2 ${srcdir}/rpe_ex.vgl2
	./tester xev_ex1.out ${srcdir}/xev_ex1.vgl
	./tester xev_ex2.out ${srcdir}/xev_ex2.vgl
	./tester xev_ex3.out ${srcdir}/xev_ex3.vgl
	./sumup
