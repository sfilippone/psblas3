include ../../Make.inc

#FCOPT= $(FCOPT) 
F90_PSDOBJS= dcsdp90.o dcssm90.o dcssm90v.o dfixcoo.o dipcoo2csr.o  dipcsr2coo.o\
	  dcsprt90.o dspgtdiag.o dspinfo.o dspgtrow.o dspscal.o imsr.o imsrx.o \
	  dsymbmm90.o dnumbmm90.o drwextd.o dtransp90.o smmp.o dcsmm90.o dcsmv90.o\
	  dcsrws90.o psdneigh.o psbdcoins.o string_impl.o dcsnmi90.o


LIBDIR= ../../lib
INCLUDES=-I$(LIBDIR) -I..
LIBNAME=$(LIBDIR)/$(F90LIB)
HERE=.

INCDIRS=-I. -I.. -I$(LIBDIR)

lib: $(F90_PSDOBJS) 
	ar -cur $(LIBNAME) $(F90_PSDOBJS)
	ranlib $(LIBNAME)

#$(F90_PSDOBJS): $(MODS) 

.f.o:
	$(F90) $(FCOPT) $(INCDIRS)  -c $<

veryclean: clean
	/bin/rm -f $(LIBNAME)

clean: 
	/bin/rm -f $(F90_PSDOBJS) $(LOCAL_MODS)
