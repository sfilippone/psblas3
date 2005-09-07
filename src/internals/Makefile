include ../../Make.inc

FOBJS = psi_compute_size.o psi_crea_bnd_elem.o psi_crea_index.o \
	psi_crea_ovr_elem.o psi_dl_check.o \
	psi_exist_ovr_elem.o psi_gthsct.o \
	psi_list_search.o psi_sort_dl.o srtlist.o
COBJS = avltree.o

MPFOBJS = psi_dswapdata.o psi_dswaptran.o psi_iswapdata.o \
	psi_iswaptran.o psi_extrct_dl.o psi_desc_index.o
INCDIRS = ../../lib .

lib: mpfobjs $(FOBJS) $(COBJS)


mpfobjs: 
	(make $(MPFOBJS) F90="$(MPF90)" FC="$(MPF90)" FCOPT="$(F90COPT)")


clean:
	/bin/rm -f $(MPFOBJS) $(FOBJS) $(COBJS)
