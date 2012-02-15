include ../Make.inc


LIBDIR=../lib

HERE=.

BASEOBJS= psb_blockpart_mod.o psb_metispart_mod.o \
      psb_hbio_mod.o psb_mmio_mod.o psb_mat_dist_mod.o \
      psb_renum_mod.o psb_gps_mod.o psb_d_genmat_mod.o
IMPLOBJS= psb_s_hbio_impl.o psb_d_hbio_impl.o \
	psb_c_hbio_impl.o psb_z_hbio_impl.o \
	psb_s_mmio_impl.o psb_d_mmio_impl.o \
	psb_c_mmio_impl.o psb_z_mmio_impl.o \
	psb_s_mat_dist_impl.o psb_d_mat_dist_impl.o \
	psb_c_mat_dist_impl.o psb_z_mat_dist_impl.o \
	psb_s_renum_impl.o psb_d_renum_impl.o \
	psb_c_renum_impl.o psb_z_renum_impl.o \
	psb_d_genmat_impl.o
MODOBJS=psb_util_mod.o $(BASEOBJS)
COBJS=psb_amd_order.o
OBJS=$(MODOBJS) $(IMPLOBJS) $(COBJS)
LOCAL_MODS=$(MODOBJS:.o=$(.mod))
LIBNAME=$(UTILLIBNAME)
FINCLUDES=$(FMFLAG). $(FMFLAG)$(LIBDIR) 


lib:  $(HERE)/$(LIBNAME) 
	/bin/cp -p $(CPUPDFLAG) $(HERE)/$(LIBNAME) $(LIBDIR)
	/bin/cp -p $(CPUPDFLAG) *$(.mod) $(LIBDIR)

$(HERE)/$(LIBNAME): $(OBJS)
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)

$(OBJS): $(LIBDIR)/$(BASEMODNAME)$(.mod)
psb_util_mod.o: $(BASEOBJS)

$(IMPLOBJS): $(BASEOBJS)

veryclean: clean
	/bin/rm -f $(HERE)/$(LIBNAME)

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

veryclean: clean
