include ../../Make.inc

MODULES = psbn_base_mat_mod.o psbn_d_base_mat_mod.o psbn_d_coo_impl.o 
#	psbn_csr_mat.o 


LIBMOD=
OBJS = 
LIBDIR=..
CINCLUDES=-I. 
FINCLUDES=$(FMFLAG)$(LIBDIR) $(FMFLAG). $(FIFLAG).


lib:  $(MODULES)  $(OBJS)  $(LIBMOD)
	$(AR) $(LIBDIR)/$(LIBNAME) $(MODULES) $(OBJS) $(MPFOBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME) #	/bin/cp -p $(LIBMOD) $(LIBDIR)
	/bin/cp -p *$(.mod) $(LIBDIR)

psbn_mat_mod.o:	psbn_base_mat_mod.o
psbn_coo_mat.o psbn_csr_mat.o:	psbn_d_base_mat_mod.o

psbn_d_mat_impl.o: psbn_d_base_mat_mod.o

clean:
	/bin/rm -f $(MODULES) $(OBJS) $(MPFOBJS) *$(.mod)

veryclean: clean

