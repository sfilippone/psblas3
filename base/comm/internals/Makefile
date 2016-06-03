include ../../../Make.inc

FOBJS = psi_iovrl_restr.o psi_iovrl_save.o psi_iovrl_upd.o \
	psi_sovrl_restr.o psi_sovrl_save.o psi_sovrl_upd.o \
	psi_dovrl_restr.o psi_dovrl_save.o psi_dovrl_upd.o \
	psi_covrl_restr.o psi_covrl_save.o psi_covrl_upd.o \
	psi_zovrl_restr.o psi_zovrl_save.o psi_zovrl_upd.o 

MPFOBJS = psi_dswapdata.o psi_dswaptran.o\
	psi_sswapdata.o  psi_sswaptran.o \
	psi_iswapdata.o psi_iswaptran.o \
	psi_cswapdata.o psi_cswaptran.o  \
	psi_zswapdata.o psi_zswaptran.o  
LIBDIR=../..
INCDIR=../..
MODDIR=../../modules
FINCLUDES=$(FMFLAG). $(FMFLAG)$(MODDIR) $(FMFLAG)$(INCDIR)
CINCLUDES=-I. 

lib: mpfobjs $(FOBJS) $(MPFOBJS) 
	$(AR) $(LIBDIR)/$(LIBNAME) $(MPFOBJS) $(MPFOBJS2) $(FOBJS) $(FOBJS2) $(COBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)

$(FOBJS) $(FBOJS2): $(MODDIR)/psi_mod.o
mpfobjs: 
	(make $(MPFOBJS) F90="$(MPF90)" FC="$(MPF90)" FCOPT="$(F90COPT)")	

clean:
	/bin/rm -f $(MPFOBJS) $(FOBJS) $(COBJS) $(FOBJS2) $(MPFOBJS2) *$(.mod)

veryclean: clean
