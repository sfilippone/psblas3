include ../../Make.inc

FOBJS = psi_crea_bnd_elem.o psi_crea_index.o \
        psi_crea_ovr_elem.o psi_bld_tmpovrl.o \
	psi_bld_tmphalo.o psi_sort_dl.o \
	psi_indx_map_fnd_owner.o \
        psi_desc_impl.o psi_hash_impl.o psi_srtlist.o \
        psi_bld_glb_dep_list.o psi_xtr_loc_dl.o

MPFOBJS = psi_desc_index.o psi_fnd_owner.o psi_a2a_fnd_owner.o \
	psi_graph_fnd_owner.o psi_adjcncy_fnd_owner.o psi_symm_dep_list.o

LIBDIR=..
INCDIR=..
MODDIR=../modules
FINCLUDES=$(FMFLAG). $(FMFLAG)$(MODDIR) $(FMFLAG)$(INCDIR)
CINCLUDES=-I. 

objs: mpfobjs $(FOBJS) $(FOBJS2) $(COBJS) $(MPFOBJS2) $(MPFOBJS)
lib: objs
	$(AR) $(LIBDIR)/$(LIBNAME) $(MPFOBJS) $(MPFOBJS2) $(FOBJS) $(FOBJS2) $(COBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)

$(FOBJS) $(FBOJS2): $(MODDIR)/psi_mod.o
mpfobjs: 
	(make $(MPFOBJS) FC="$(MPFC)")	

clean:
	/bin/rm -f $(MPFOBJS) $(FOBJS) $(COBJS) $(FOBJS2) $(MPFOBJS2) *$(.mod)

veryclean: clean
