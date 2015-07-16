include ../../Make.inc

FOBJS = psi_compute_size.o psi_crea_bnd_elem.o psi_crea_index.o \
        psi_crea_ovr_elem.o psi_bld_tmpovrl.o psi_dl_check.o \
	psi_sort_dl.o \
	psi_bld_tmphalo.o\
	psi_sort_dl.o \
        psi_desc_impl.o \
	psi_iovrl_restr.o psi_iovrl_save.o psi_iovrl_upd.o \
	psi_sovrl_restr.o psi_sovrl_save.o psi_sovrl_upd.o \
	psi_dovrl_restr.o psi_dovrl_save.o psi_dovrl_upd.o \
	psi_covrl_restr.o psi_covrl_save.o psi_covrl_upd.o \
	psi_zovrl_restr.o psi_zovrl_save.o psi_zovrl_upd.o 

FOBJS2 = psi_exist_ovr_elem.o psi_list_search.o srtlist.o 
#COBJS = avltree.o srcht.o

MPFOBJS = psi_dswapdata.o psi_dswaptran.o\
	psi_sswapdata.o  psi_sswaptran.o \
	psi_iswapdata.o psi_iswaptran.o \
	psi_cswapdata.o psi_cswaptran.o  \
	psi_zswapdata.o psi_zswaptran.o  \
	psi_desc_index.o psi_extrct_dl.o \
	psi_fnd_owner.o psb_indx_map_fnd_owner.o
LIBDIR=..
INCDIR=..
MODDIR=../modules
FINCLUDES=$(FMFLAG). $(FMFLAG)$(MODDIR) $(FMFLAG)$(INCDIR)
CINCLUDES=-I. 

lib: mpfobjs $(FOBJS) $(FOBJS2) $(COBJS) $(MPFOBJS2) $(MPFOBJS) 
	$(AR) $(LIBDIR)/$(LIBNAME) $(MPFOBJS) $(MPFOBJS2) $(FOBJS) $(FOBJS2) $(COBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)

$(FOBJS) $(FBOJS2): $(MODDIR)/psi_mod.o
mpfobjs: 
	(make $(MPFOBJS) F90="$(MPF90)" FC="$(MPF90)" FCOPT="$(F90COPT)")	
	(make $(FOBJS2) F90="$(MPF77)" FC="$(MPF77)" FCOPT="$(FCOPT)")
clean:
	/bin/rm -f $(MPFOBJS) $(FOBJS) $(COBJS) $(FOBJS2) $(MPFOBJS2) *$(.mod)

veryclean: clean
