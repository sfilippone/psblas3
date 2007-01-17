include ../Make.inc


HERE=.
LIBDIR=../lib

OBJS=psb_krylov_mod.o psb_dcgstab.o psb_dcg.o psb_dcgs.o \
	psb_dbicg.o psb_dcgstabl.o psb_dgmresr.o\
	psb_zcgstab.o psb_zcgs.o


LIBMOD=psb_krylov_mod$(.mod)
LOCAL_MODS=$(LIBMOD)
LIBNAME=$(METHDLIBNAME)	

INCDIRS=-I. -I$(LIBDIR) 

lib:  $(OBJS) 
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)
	/bin/cp -p $(HERE)/$(LIBNAME) $(LIBDIR)
	/bin/cp -p $(LIBMOD) $(LIBDIR)

veryclean: clean
	/bin/rm -f $(HERE)/$(LIBNAME)

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

veryclean: clean
