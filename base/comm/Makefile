include ../../Make.inc

OBJS =  psb_dgather.o psb_dhalo.o psb_dovrl.o \
	psb_sgather.o psb_shalo.o psb_sovrl.o \
	psb_igather.o psb_ihalo.o psb_iovrl.o \
	psb_lgather.o psb_lhalo.o psb_lovrl.o \
	psb_cgather.o psb_chalo.o psb_covrl.o \
	psb_zgather.o psb_zhalo.o psb_zovrl.o \
	psb_dgather_a.o psb_dhalo_a.o psb_dovrl_a.o \
	psb_sgather_a.o psb_shalo_a.o psb_sovrl_a.o \
	psb_mgather_a.o psb_mhalo_a.o psb_movrl_a.o \
	psb_egather_a.o psb_ehalo_a.o psb_eovrl_a.o \
	psb_cgather_a.o psb_chalo_a.o psb_covrl_a.o \
	psb_zgather_a.o psb_zhalo_a.o psb_zovrl_a.o 


MPFOBJS=psb_dscatter.o psb_zscatter.o \
	psb_iscatter.o psb_lscatter.o \
	psb_cscatter.o psb_sscatter.o \
	psb_dscatter_a.o psb_zscatter_a.o \
	psb_mscatter_a.o psb_escatter_a.o \
	psb_cscatter_a.o psb_sscatter_a.o \
	psb_dspgather.o psb_sspgather.o \
	psb_zspgather.o psb_cspgather.o
LIBDIR=..
INCDIR=..
MODDIR=../modules
FINCLUDES=$(FMFLAG). $(FMFLAG)$(MODDIR) $(FMFLAG)$(INCDIR) 

objs: interns mpfobjs $(OBJS)

lib: objs
	$(MAKE) -C internals lib LIBNAME=$(LIBNAME)
	$(AR) $(LIBDIR)/$(LIBNAME) $(MPFOBJS) $(OBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)

interns:
	$(MAKE) -C internals objs
mpfobjs: 
	$(MAKE) $(MPFOBJS) FC="$(MPFC)" 


clean:
	cd internals && $(MAKE) clean
	/bin/rm -f $(MPFOBJS) $(OBJS)

veryclean: clean
