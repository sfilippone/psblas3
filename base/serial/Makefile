include ../../Make.inc


FOBJS = psb_lsame.o psi_serial_impl.o psb_sort_impl.o \
	psb_srwextd.o psb_drwextd.o psb_crwextd.o psb_zrwextd.o \
	psb_ssymbmm.o psb_dsymbmm.o psb_csymbmm.o psb_zsymbmm.o \
	psb_snumbmm.o psb_dnumbmm.o psb_cnumbmm.o psb_znumbmm.o \
	psb_sgeprt.o psb_dgeprt.o psb_cgeprt.o psb_zgeprt.o\
	psb_spdot_srtd.o psb_aspxpby.o psb_spge_dot.o 

LIBDIR=..
MODDIR=../modules
FINCLUDES=$(FMFLAG)$(LIBDIR) $(FMFLAG)$(MODDIR) $(FMFLAG).

lib: auxd  f77d impld lib1 $(FOBJS)
	$(AR) $(LIBDIR)/$(LIBNAME) $(FOBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)

#cood csrd jadd

lib1:	$(FOBJS)


auxd: 
	(cd aux; $(MAKE) lib)

f77d:
	(cd f77; $(MAKE) lib)
impld:
	(cd impl; $(MAKE) lib)

clean:
	/bin/rm -f $(FOBJS) *$(.mod)
	(cd aux; $(MAKE) clean)
	(cd f77; $(MAKE) clean)
	(cd impl; $(MAKE) clean)

veryclean: clean 
