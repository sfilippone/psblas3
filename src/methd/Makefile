include ../../Make.inc


LIBDIR=../../lib

HERE=.
OBJS= psb_dcgstab.o psb_dcg.o psb_dcgs.o \
	psb_dbicg.o psb_dcgstabl.o psb_dgmresr.o

INCDIRS=-I. -I.. -I$(LIBDIR)

lib:  $(OBJS) 
	$(AR) $(LIBDIR)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)

#$(F90OBJS): $(MODS)  


veryclean: clean
	/bin/rm -f $(LIBNAME)

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

veryclean: clean
