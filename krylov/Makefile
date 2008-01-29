include ../Make.inc


HERE=.
LIBDIR=../lib

MODOBJS= psb_krylov_mod.o
F90OBJS= psb_dcgstab.o psb_dcg.o psb_dcgs.o \
	psb_dbicg.o psb_dcgstabl.o psb_drgmres.o\
	psb_zcgstab.o psb_zcg.o psb_zcgs.o  \
	psb_zbicg.o psb_zcgstabl.o psb_zrgmres.o
OBJS=$(F90OBJS) $(MODOBJS)

LIBMOD=psb_krylov_mod$(.mod)
LOCAL_MODS=$(LIBMOD)
LIBNAME=$(METHDLIBNAME)	

FINCLUDES=$(FMFLAG)$(LIBDIR) $(FMFLAG).

lib:  $(OBJS) 
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)
	/bin/cp -p $(HERE)/$(LIBNAME) $(LIBDIR)
	/bin/cp -p $(LIBMOD) $(LIBDIR)

$(F90OBJS): $(MODOBJS)
$(OBJS): $(LIBDIR)/psb_prec_mod$(.mod)  $(LIBDIR)/psb_base_mod$(.mod)
veryclean: clean
	/bin/rm -f $(HERE)/$(LIBNAME)

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

veryclean: clean
