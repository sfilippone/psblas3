include ../../Make.inc

#FCOPT=-O2
F90_PSDOBJS= psb_ddot.o psb_damax.o psb_dasum.o psb_daxpby.o\
	 psb_dnrm2.o psb_dnrmi.o psb_dspmm.o psb_dspsm.o\

LIBDIR=../../lib
HERE=.
LIBNAME=$(LIBDIR)/$(F90LIB)


INCDIRS=-I. -I.. -I$(LIBDIR)


lib: $(F90_PSDOBJS) 
	(cd INTERNALS; make lib LIBDIR=../$(LIBDIR) LIBNAME=$(LIBNAME)) 
	ar -cur $(LIBNAME) $(F90_PSDOBJS) 
	ranlib $(LIBNAME)

#$(F90_PSDOBJS): $(MODS)  

veryclean: clean
	/bin/rm -f $(LIBNAME)

clean: 
	/bin/rm -f $(F90_PSDOBJS) $(LOCAL_MODS)
