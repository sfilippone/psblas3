include ../../Make.inc

#FCOPT=-O2
OBJS= psb_ddot.o psb_damax.o psb_dasum.o psb_daxpby.o\
	 psb_dnrm2.o psb_dnrmi.o psb_dspmm.o psb_dspsm.o\
	pdtreecomb.o

LIBDIR=../../lib
HERE=.

INCDIRS=-I. -I.. -I$(LIBDIR)


lib: $(OBJS) 
	$(AR) $(LIBDIR)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)

#$(F90_PSDOBJS): $(MODS)  

veryclean: clean
	/bin/rm -f $(LIBNAME)

clean: 
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

veryclean: clean
