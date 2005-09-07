include ../../Make.inc


LIBDIR=../../lib/
LIBNAME=$(LIBDIR)/$(F90LIB)

HERE=.
F90OBJS= f90_dcgstab.o f90_dcg.o f90_dcgs.o \
	f90_dbicg.o f90_dcgstabl.o f90_zcgstab.o f90_dgmresr.o

INCDIRS=-I. -I.. -I$(LIBDIR)

lib:  $(F90OBJS) 
	ar -cur $(LIBNAME)  $(F90OBJS) 
	ranlib $(LIBNAME)

#$(F90OBJS): $(MODS)  


veryclean: clean
	/bin/rm -f $(LIBNAME)

clean:
	/bin/rm -f $(F90OBJS) $(LOCAL_MODS)
