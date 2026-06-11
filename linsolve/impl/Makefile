include ../../Make.inc

LIBDIR=../../lib
INCDIR=../../include
MODDIR=../../modules
HERE=..
OBJS=psb_dkrylov.o psb_skrylov.o psb_ckrylov.o psb_zkrylov.o \
	psb_drichardson.o psb_srichardson.o psb_crichardson.o psb_zrichardson.o \
        psb_dcgstab.o psb_dcg.o psb_dfcg.o psb_dgcr.o psb_dcgs.o \
	psb_dbicg.o psb_dcgstabl.o psb_drgmres.o\
	psb_scgstab.o psb_scg.o psb_sfcg.o psb_sgcr.o psb_scgs.o \
	psb_sbicg.o psb_scgstabl.o psb_srgmres.o\
	psb_ccgstab.o psb_ccg.o psb_cfcg.o psb_cgcr.o psb_ccgs.o  \
	psb_cbicg.o psb_ccgstabl.o psb_crgmres.o\
	psb_zcgstab.o psb_zcg.o psb_zfcg.o psb_zgcr.o psb_zcgs.o  \
	psb_zbicg.o psb_zcgstabl.o psb_zrgmres.o \
	psb_dminres.o psb_sminres.o psb_cminres.o psb_zminres.o 

LIBNAME=$(METHDLIBNAME)
COBJS=
FINCLUDES=$(FMFLAG).. $(FMFLAG)$(MODDIR)

objs:  $(OBJS)
lib:  objs
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)

veryclean: clean

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)
