include ../../Make.inc

FOBJS = psi_compute_size.o psi_crea_bnd_elem.o psi_crea_index.o \
        psi_crea_ovr_elem.o psi_bld_tmpovrl.o psi_dl_check.o \
	psi_bld_tmphalo.o\
	psi_sort_dl.o \
        psi_desc_impl.o

FOBJS2 = psi_exist_ovr_elem.o psi_list_search.o srtlist.o 

MPFOBJS = psi_desc_index.o psi_extrct_dl.o \
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
