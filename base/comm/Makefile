include ../../Make.inc

OBJS = psb_dgather.o psb_dhalo.o psb_dovrl.o \
	psb_ihalo.o psb_zgather.o psb_zhalo.o psb_zovrl.o 

MPFOBJS = psb_dscatter.o  psb_zscatter.o
LIBDIR = ..
MODDIR = ../modules
INCDIRS = -I $(LIBDIR) -I $(MODDIR) -I . 

lib: mpfobjs $(OBJS)
	$(AR) $(LIBDIR)/$(LIBNAME) $(MPFOBJS) $(OBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)


mpfobjs: 
	(make $(MPFOBJS) F90="$(MPF90)" FC="$(MPF90)" FCOPT="$(F90COPT)")


clean:
	/bin/rm -f $(MPFOBJS) $(OBJS)

veryclean: clean
