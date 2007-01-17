include ../Make.inc


LIBDIR=../lib

HERE=.

BASEOBJS= psb_blockpart_mod.o psb_metispart_mod.o \
      psb_hbio_mod.o psb_mmio_mod.o psb_mat_dist_mod.o \
      psb_read_mat_mod.o 
MODOBJ=psb_util_mod.o
OBJS=$(BASEOBJS) $(MODOBJ)
LIBMOD=psb_util_mod$(.mod)
LOCAL_MODS=$(OBJS:.o=$(.mod))
LIBNAME=$(UTILLIBNAME)
INCDIRS=-I.  -I$(LIBDIR)

lib:  $(OBJS) 
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)
	/bin/cp -p $(HERE)/$(LIBNAME) $(LIBDIR)
	/bin/cp -p $(LIBMOD) $(LIBDIR)


psb_util_mod.o: $(BASEOBJ)

veryclean: clean
	/bin/rm -f $(HERE)/$(LIBNAME)

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

veryclean: clean
