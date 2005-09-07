include ../../Make.inc


LIBDIR=../../lib/
LIBNAME=$(LIBDIR)/$(F90LIB)

HERE=.
MPFOBJS=dcslu.o psbdbldaggrmat.o 
F90OBJS= dcsrsetup.o dcsrlu.o f90_psdprec.o \
	dprecbld.o zprecbld.o gps.o psdprecfree.o dprecset.o \
	 psbdgenaggrmap.o  $(MPFOBJS) 
#dcoocp.o dcoocpadd.o dcoofact.o dcoolu.o dcooluadd.o\

COBJS=fort_slu_impl.o
INCDIRS=-I. -I.. -I$(LIBDIR) 

OBJS=$(F90OBJS) $(COBJS)

lib:  mpobjs $(OBJS) 
	ar -cur $(LIBNAME)  $(OBJS) 
	ranlib $(LIBNAME)

#$(F90OBJS): $(MODS)  
mpobjs: 
	(make $(MPFOBJS) F90="$(MPF90)" F90COPT="$(F90COPT)")

veryclean: clean
	/bin/rm -f $(LIBNAME)

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)
