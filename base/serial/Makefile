include ../../Make.inc


FOBJS = psb_lsame.o psi_serial_impl.o psi_impl.o psb_sort_impl.o \
	psb_ssymbmm.o psb_dsymbmm.o psb_csymbmm.o psb_zsymbmm.o \
	psb_snumbmm.o psb_dnumbmm.o psb_cnumbmm.o psb_znumbmm.o \
	psb_srwextd.o psb_drwextd.o psb_crwextd.o psb_zrwextd.o

LIBDIR=..
MODDIR=../modules
FINCLUDES=$(FMFLAG)$(LIBDIR) $(FMFLAG)$(MODDIR) $(FMFLAG).

lib: auxd  f77d f03d lib1
	$(AR) $(LIBDIR)/$(LIBNAME) $(FOBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)

#cood csrd jadd

lib1:	$(FOBJS)


auxd: 
	(cd aux; make lib)

f77d:
	(cd f77; make lib)
f03d:
	(cd f03; make lib)

clean:
	/bin/rm -f $(FOBJS) *$(.mod)
	(cd aux; make clean)
	(cd f77; make clean)
	(cd f03; make clean)

veryclean: clean 
