include ../../Make.inc

#FCOPT=-O2
OBJS= psb_ddot.o psb_damax.o psb_dasum.o psb_daxpby.o\
	 psb_dnrm2.o psb_dnrmi.o psb_dspmm.o psb_dspsm.o\
	psb_dspnrm1.o \
	psb_zamax.o psb_zasum.o psb_zaxpby.o psb_zdot.o \
	psb_znrm2.o psb_znrmi.o psb_zspmm.o psb_zspsm.o\
	psb_saxpby.o psb_sdot.o  psb_sasum.o psb_samax.o\
	psb_snrm2.o psb_snrmi.o  psb_sspmm.o psb_sspsm.o\
	psb_camax.o psb_casum.o psb_caxpby.o psb_cdot.o \
	psb_cnrm2.o psb_cnrmi.o psb_cspmm.o psb_cspsm.o


LIBDIR=..
MODDIR=../modules
FINCLUDES=$(FMFLAG)$(LIBDIR) $(FMFLAG)$(MODDIR) $(FMFLAG).


lib: $(OBJS) 
	$(AR) $(LIBDIR)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)

#$(F90_PSDOBJS): $(MODS)  

veryclean: clean
	/bin/rm -f $(LIBNAME)

clean: 
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

veryclean: clean
