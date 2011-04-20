include ../Make.inc
#
# Libraries used
#
LIBDIR=../lib/
PSBLAS_LIB= -L$(LIBDIR) -lpsb_base
#-lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base 
LDLIBS=$(PSBLDLIBS)
#
# Compilers and such
#
CCOPT= -g
FINCLUDES=$(FMFLAG)$(LIBDIR) $(FMFLAG).
HERE=.
LIBMOD=*$(.mod)

LIBDIR=../lib

EXEDIR=./runs

OBJS=psb_d_ell_impl.o psb_d_ell_mat_mod.o \
	rsb_z_mod.o psb_z_rsb_mat_mod.o \
	rsb_c_mod.o psb_c_rsb_mat_mod.o \
	rsb_s_mod.o psb_s_rsb_mat_mod.o \
       	rsb_d_mod.o psb_d_rsb_mat_mod.o

LIBNAME=libpsb_opt.a

lib: $(OBJS)
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)
	/bin/cp -p $(CPUPDFLAG) $(HERE)/$(LIBNAME) $(LIBDIR)
	/bin/cp -p $(CPUPDFLAG)  *$(.mod) $(LIBDIR)


libpsb_opt.a: $(OBJS)
	ar cur libpsb_opt.a $(OBJS)

psb_d_ell_impl.o: psb_d_ell_mat_mod.o
psb_d_rsb_mat_mod.o:  rsb_d_mod.o
psb_z_rsb_mat_mod.o:  rsb_z_mod.o
psb_c_rsb_mat_mod.o:  rsb_c_mod.o
psb_s_rsb_mat_mod.o:  rsb_s_mod.o

RSBD2Z=sed 's/rsb_d_/rsb_z_/g;s/psb_d_/psb_z_/g;s/real(psb_dpk_)/complex(psb_dpk_)/g;s/real(c_double)/complex(c_double)/g;s/complex(psb_dpk_)\(.*\)csnmi_res/real(psb_dpk_)\1csnmi_res/g' 
RSBD2S=sed 's/rsb_d_/rsb_s_/g;s/psb_d_/psb_s_/g;s/real(psb_dpk_)/real(psb_spk_)/g;s/real(c_double)/real(c_float)/g;s/real(psb_dpk_)\(.*\)csnmi_res/real(psb_spk_)\1csnmi_res/g' 
RSBD2C=sed 's/rsb_d_/rsb_c_/g;s/psb_d_/psb_c_/g;s/real(psb_dpk_)/complex(psb_spk_)/g;s/real(c_double)/complex(c_float)/g;s/complex(psb_spk_)\(.*\)csnmi_res/real(psb_spk_)\1csnmi_res/g' 

psb_z_rsb_mat_mod.F90: psb_d_rsb_mat_mod.F90 Makefile
	$(RSBD2Z) $< > $@

rsb_z_mod.f90: rsb_d_mod.f90 Makefile
	$(RSBD2Z) $< > $@

psb_c_rsb_mat_mod.F90: psb_d_rsb_mat_mod.F90 Makefile
	$(RSBD2C) $< > $@

rsb_c_mod.f90: rsb_d_mod.f90 Makefile
	$(RSBD2C) $< > $@

psb_s_rsb_mat_mod.F90: psb_d_rsb_mat_mod.F90 Makefile
	$(RSBD2S) $< > $@

rsb_s_mod.f90: rsb_d_mod.f90 Makefile
	$(RSBD2S) $< > $@


clean: 
	/bin/rm -f  $(OBJS) *$(.mod) 

verycleanlib: 
	(cd ../..; make veryclean)




