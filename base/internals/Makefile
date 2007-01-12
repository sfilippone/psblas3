include ../../Make.inc

FOBJS = psi_compute_size.o psi_crea_bnd_elem.o psi_crea_index.o \
        psi_crea_ovr_elem.o psi_dl_check.o \
	psi_gthsct.o \
	psi_sort_dl.o \
	psi_gthsct.o psi_ldsc_pre_halo.o\
	psi_sort_dl.o psi_idx_cnv.o psi_idx_ins_cnv.o psi_fnd_owner.o
FOBJS2 = psi_exist_ovr_elem.o psi_list_search.o srtlist.o 
COBJS = avltree.o srcht.o

MPFOBJS = psi_dswapdata.o psi_dswaptran.o psi_iswapdata.o \
	psi_iswaptran.o psi_desc_index.o \
 	psi_zswapdata.o psi_zswaptran.o 
MPFOBJS2 = psi_extrct_dl.o
LIBDIR = ..
MODDIR = ../modules
INCDIRS = -I $(LIBDIR) -I $(MODDIR) -I . 

lib: mpfobjs $(FOBJS) $(FOBJS2) $(COBJS) $(MPFOBJS2)
	$(AR) $(LIBDIR)/$(LIBNAME) $(MPFOBJS) $(MPFOBJS2) $(FOBJS) $(FOBJS2) \
             $(COBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)


mpfobjs: psi_gthsct.o
	(make $(MPFOBJS) F90="$(MPF90)" FC="$(MPF90)" FCOPT="$(F90COPT)")	
	(make $(MPFOBJS2) F90="$(MPF77)" FC="$(MPF77)" FCOPT="$(FCOPT)")
	(make $(FOBJS2) F90="$(MPF77)" FC="$(MPF77)" FCOPT="$(FCOPT)")
clean:
	/bin/rm -f $(MPFOBJS) $(FOBJS) $(COBJS) $(FOBJS2) $(MPFOBJS2)

veryclean: clean
