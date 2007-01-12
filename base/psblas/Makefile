include ../../Make.inc

#FCOPT=-O2
OBJS= psb_ddot.o psb_damax.o psb_dasum.o psb_daxpby.o\
	 psb_dnrm2.o psb_dnrmi.o psb_dspmm.o psb_dspsm.o\
	pdtreecomb.o \
	psb_zamax.o psb_zasum.o psb_zaxpby.o psb_zdot.o \
	psb_znrm2.o psb_znrmi.o psb_zspmm.o psb_zspsm.o

LIBDIR = ..
MODDIR = ../modules
INCDIRS = -I $(LIBDIR) -I $(MODDIR) -I . 


lib: $(OBJS) 
	$(AR) $(LIBDIR)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)

#$(F90_PSDOBJS): $(MODS)  

veryclean: clean
	/bin/rm -f $(LIBNAME)

clean: 
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

veryclean: clean
