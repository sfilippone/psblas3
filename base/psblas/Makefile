include ../../Make.inc

#FCOPT=-O2
OBJS= psb_ddot.o psb_damax.o psb_dasum.o psb_daxpby.o psb_dqrfact.o\
    psb_dprod.o psb_dnrm2.o psb_dnrmi.o psb_dspmm.o psb_dspsm.o\
	psb_sspnrm1.o psb_dspnrm1.o psb_cspnrm1.o psb_zspnrm1.o \
	psb_zamax.o psb_zasum.o psb_zaxpby.o psb_zdot.o \
	psb_znrm2.o psb_znrmi.o psb_zspmm.o psb_zspsm.o\
	psb_saxpby.o psb_sdot.o  psb_sasum.o psb_samax.o\
	psb_snrm2.o psb_snrmi.o  psb_sspmm.o psb_sspsm.o\
	psb_camax.o psb_casum.o psb_caxpby.o psb_cdot.o \
	psb_cnrm2.o psb_cnrmi.o psb_cspmm.o psb_cspsm.o \
	psb_cmlt_vect.o psb_dmlt_vect.o psb_zmlt_vect.o psb_smlt_vect.o\
	psb_cdiv_vect.o psb_ddiv_vect.o psb_zdiv_vect.o psb_sdiv_vect.o\
	psb_cinv_vect.o psb_dinv_vect.o psb_zinv_vect.o psb_sinv_vect.o\
	psb_dcmp_vect.o psb_scmp_vect.o psb_ccmp_vect.o psb_zcmp_vect.o\
	psb_cabs_vect.o psb_dabs_vect.o psb_sabs_vect.o \
	psb_zabs_vect.o psb_cgetmatinfo.o psb_dgetmatinfo.o psb_sgetmatinfo.o \
	psb_zgetmatinfo.o

LIBDIR=..
INCDIR=..
MODDIR=../modules
FINCLUDES=$(FMFLAG). $(FMFLAG)$(MODDIR) $(FMFLAG)$(INCDIR)


objs: $(OBJS)

lib: objs
	$(AR) $(LIBDIR)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)

#$(F90_PSDOBJS): $(MODS)

veryclean: clean
	/bin/rm -f $(LIBNAME)

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

veryclean: clean
