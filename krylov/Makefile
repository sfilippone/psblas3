include ../Make.inc


HERE=.
LIBDIR=../lib
INCDIR=../include

MODOBJS= psb_base_krylov_conv_mod.o \
	psb_s_krylov_conv_mod.o psb_c_krylov_conv_mod.o \
	psb_d_krylov_conv_mod.o psb_z_krylov_conv_mod.o \
	psb_krylov_mod.o
F90OBJS=psb_dkrylov.o psb_skrylov.o psb_ckrylov.o psb_zkrylov.o \
        psb_dcgstab.o psb_dcg.o psb_dfcg.o psb_dcgr.o psb_dcgs.o \
	psb_dbicg.o psb_dcgstabl.o psb_drgmres.o\
	psb_scgstab.o psb_scg.o psb_sfcg.o psb_scgr.o psb_scgs.o \
	psb_sbicg.o psb_scgstabl.o psb_srgmres.o\
	psb_ccgstab.o psb_ccg.o psb_cfcg.o psb_ccgr.o psb_ccgs.o  \
	psb_cbicg.o psb_ccgstabl.o psb_crgmres.o\
	psb_zcgstab.o psb_zcg.o psb_zfcg.o psb_zcgr.o psb_zcgs.o  \
	psb_zbicg.o psb_zcgstabl.o psb_zrgmres.o
OBJS=$(F90OBJS) $(MODOBJS)

LOCAL_MODS=$(MODOBJS:.o=$(.mod))
LIBNAME=$(METHDLIBNAME)	

FINCLUDES=$(FMFLAG). $(FMFLAG)$(INCDIR) 

lib: $(OBJS) 
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)
	/bin/cp -p $(CPUPDFLAG) $(HERE)/$(LIBNAME) $(LIBDIR)
	/bin/cp -p $(CPUPDFLAG) *$(.mod) $(INCDIR)


psb_s_krylov_conv_mod.o psb_c_krylov_conv_mod.o psb_d_krylov_conv_mod.o psb_z_krylov_conv_mod.o:  psb_base_krylov_conv_mod.o
psb_krylov_conv_mod.o: psb_s_krylov_conv_mod.o psb_c_krylov_conv_mod.o psb_d_krylov_conv_mod.o psb_z_krylov_conv_mod.o 
$(F90OBJS): $(MODOBJS)
$(OBJS): $(INCDIR)/$(PRECMODNAME)$(.mod)  $(INCDIR)/$(BASEMODNAME)$(.mod)

veryclean: clean
	/bin/rm -f $(HERE)/$(LIBNAME)

clean:
	/bin/rm -f $(OBJS) *$(.mod)

