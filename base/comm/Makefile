include ../../Make.inc

OBJS =  psb_dgather.o psb_dhalo.o psb_dovrl.o \
	psb_sgather.o psb_shalo.o psb_sovrl.o \
	psb_igather.o psb_ihalo.o psb_iovrl.o \
	psb_cgather.o psb_chalo.o psb_covrl.o \
	psb_zgather.o psb_zhalo.o psb_zovrl.o 


MPFOBJS=psb_dscatter.o psb_zscatter.o psb_iscatter.o psb_cscatter.o psb_sscatter.o\
	psb_dspgather.o psb_sspgather.o psb_zspgather.o psb_cspgather.o
LIBDIR=..
MODDIR=../modules
FINCLUDES=$(FMFLAG)$(LIBDIR) $(FMFLAG)$(MODDIR) $(FMFLAG).

lib: mpfobjs $(OBJS) 
	$(AR) $(LIBDIR)/$(LIBNAME) $(MPFOBJS) $(OBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)


mpfobjs: 
	(make $(MPFOBJS) F90="$(MPF90)" FC="$(MPF90)" FCOPT="$(F90COPT)")


clean:
	/bin/rm -f $(MPFOBJS) $(OBJS)

veryclean: clean
