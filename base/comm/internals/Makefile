include ../../../Make.inc

FOBJS = psi_iovrl_restr.o psi_iovrl_save.o psi_iovrl_upd.o \
	psi_lovrl_restr.o psi_lovrl_save.o psi_lovrl_upd.o \
	psi_sovrl_restr.o psi_sovrl_save.o psi_sovrl_upd.o \
	psi_dovrl_restr.o psi_dovrl_save.o psi_dovrl_upd.o \
	psi_covrl_restr.o psi_covrl_save.o psi_covrl_upd.o \
	psi_zovrl_restr.o psi_zovrl_save.o psi_zovrl_upd.o \
	psi_movrl_restr_a.o psi_movrl_save_a.o psi_movrl_upd_a.o \
	psi_eovrl_restr_a.o psi_eovrl_save_a.o psi_eovrl_upd_a.o \
	psi_sovrl_restr_a.o psi_sovrl_save_a.o psi_sovrl_upd_a.o \
	psi_dovrl_restr_a.o psi_dovrl_save_a.o psi_dovrl_upd_a.o \
	psi_covrl_restr_a.o psi_covrl_save_a.o psi_covrl_upd_a.o \
	psi_zovrl_restr_a.o psi_zovrl_save_a.o psi_zovrl_upd_a.o 

MPFOBJS = psi_dswapdata.o psi_dswaptran.o\
	psi_sswapdata.o  psi_sswaptran.o \
	psi_iswapdata.o psi_iswaptran.o \
	psi_lswapdata.o psi_lswaptran.o \
	psi_cswapdata.o psi_cswaptran.o  \
	psi_zswapdata.o psi_zswaptran.o  \
	psi_dswapdata_a.o psi_dswaptran_a.o \
	psi_sswapdata_a.o  psi_sswaptran_a.o \
	psi_mswapdata_a.o psi_mswaptran_a.o \
	psi_eswapdata_a.o psi_eswaptran_a.o \
	psi_cswapdata_a.o psi_cswaptran_a.o  \
	psi_zswapdata_a.o psi_zswaptran_a.o  
LIBDIR=../..
INCDIR=../..
MODDIR=../../modules
FINCLUDES=$(FMFLAG). $(FMFLAG)$(MODDIR) $(FMFLAG)$(INCDIR)
CINCLUDES=-I. 

objs: mpfobjs $(FOBJS) $(MPFOBJS)
lib: objs
	$(AR) $(LIBDIR)/$(LIBNAME) $(MPFOBJS) $(MPFOBJS2) $(FOBJS) $(FOBJS2) $(COBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)

$(FOBJS) $(FBOJS2): $(MODDIR)/psi_mod.o
mpfobjs: 
	(make $(MPFOBJS) FC="$(MPFC)" )	

clean:
	/bin/rm -f $(MPFOBJS) $(FOBJS) $(COBJS) $(FOBJS2) $(MPFOBJS2) *$(.mod)

veryclean: clean
