include ../Make.inc


LIBDIR=../lib

HERE=.

BASEOBJS= psb_blockpart_mod.o psb_metispart_mod.o \
      psb_hbio_mod.o psb_mmio_mod.o psb_mat_dist_mod.o \
      psb_renum_mod.o psb_gps_mod.o
IMPLOBJS= psb_hbio_impl.o psb_mmio_impl.o psb_mat_dist_impl.o psb_renum_impl.o
MODOBJS=psb_util_mod.o $(BASEOBJS)
OBJS=$(MODOBJS) $(IMPLOBJS)
LIBMOD=psb_util_mod$(.mod)
LOCAL_MODS=$(MODOBJS:.o=$(.mod))
LIBNAME=$(UTILLIBNAME)
FINCLUDES=$(FMFLAG)$(LIBDIR) $(FMFLAG).


lib:  $(LIBDIR)/$(LIBNAME) $(LIBDIR)/$(LIBMOD)

$(LIBDIR)/$(LIBNAME): $(HERE)/$(LIBNAME)
	/bin/cp -p $(CPUPDFLAG) $(HERE)/$(LIBNAME) $(LIBDIR)

$(LIBDIR)/$(LIBMOD): $(LIBMOD)
	/bin/cp -p $(CPUPDFLAG) $(LIBMOD) $(LOCAL_MODS) $(LIBDIR)

$(HERE)/$(LIBNAME): $(OBJS)
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)

psb_util_mod.o: $(BASEOBJS)

veryclean: clean
	/bin/rm -f $(HERE)/$(LIBNAME)

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

veryclean: clean
