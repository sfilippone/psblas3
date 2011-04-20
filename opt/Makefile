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

RSBD2Z=sed 's/rsb_d_/rsb_z_/g;s/psb_d_/psb_z_/g;s/real(psb_dpk_)/complex(psb_dpk_)/g;s/real(c_double)/complex(c_double)/g;s/complex(psb_dpk_)\(.*\)csnmi_res/real(psb_dpk_)\1csnmi_res/g' 

psb_z_rsb_mat_mod.F90: psb_d_rsb_mat_mod.F90 Makefile
	$(RSBD2Z) $< > $@

rsb_z_mod.f90: rsb_d_mod.f90 Makefile
	$(RSBD2Z) $< > $@


clean: 
	/bin/rm -f  $(OBJS) *$(.mod) 

verycleanlib: 
	(cd ../..; make veryclean)




