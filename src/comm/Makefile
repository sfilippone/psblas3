include ../../Make.inc

OBJS = psb_dgather.o psb_dhalo.o psb_dovrl.o \
	psb_ihalo.o

MPFOBJS = psb_dscatter.o
INCDIRS = ../../lib

lib: mpfobjs $(OBJS)


mpfobjs: 
	(make $(MPFOBJS) F90="$(MPF90)" FC="$(MPF90)" FCOPT="$(F90COPT)")


clean:
	/bin/rm -f $(MPFOBJS) $(OBJS)
