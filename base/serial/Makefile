include ../../Make.inc


FOBJS = psb_lsame.o psi_i_serial_impl.o \
	psi_s_serial_impl.o psi_d_serial_impl.o \
	psi_c_serial_impl.o psi_z_serial_impl.o \
	psb_srwextd.o psb_drwextd.o psb_crwextd.o psb_zrwextd.o \
	psb_sspspmm.o psb_dspspmm.o psb_cspspmm.o psb_zspspmm.o \
	psb_ssymbmm.o psb_dsymbmm.o psb_csymbmm.o psb_zsymbmm.o \
	psb_snumbmm.o psb_dnumbmm.o psb_cnumbmm.o psb_znumbmm.o \
	psb_sgeprt.o psb_dgeprt.o psb_cgeprt.o psb_zgeprt.o\
	psb_spdot_srtd.o psb_aspxpby.o psb_spge_dot.o\
	psb_sgelp.o psb_dgelp.o psb_cgelp.o psb_zgelp.o \
	psb_samax_s.o psb_damax_s.o psb_camax_s.o psb_zamax_s.o \
	psb_sasum_s.o psb_dasum_s.o psb_casum_s.o psb_zasum_s.o 
# psb_sort_impl.o

LIBDIR=..
INCDIR=..
MODDIR=../modules
FINCLUDES=$(FMFLAG). $(FMFLAG)$(MODDIR) $(FMFLAG)$(INCDIR) 

lib:   f77d impld sortd lib1 $(FOBJS)
	$(AR) $(LIBDIR)/$(LIBNAME) $(FOBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)

#cood csrd jadd

lib1:	$(FOBJS)


f77d:
	cd f77 && $(MAKE) lib
impld:
	cd impl && $(MAKE) lib

sortd:
	cd sort && $(MAKE) lib

clean:
	/bin/rm -f $(FOBJS) *$(.mod)
	(cd f77; $(MAKE) clean)
	(cd impl; $(MAKE) clean)
	(cd sort; $(MAKE) clean)

veryclean: clean 
