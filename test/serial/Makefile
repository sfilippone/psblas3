#
# Libraries used
#
INSTALLDIR=../..
INCDIR=$(INSTALLDIR)/include/
MODDIR=$(INSTALLDIR)/modules/
include $(INCDIR)/Make.inc.psblas
LIBDIR=$(INSTALLDIR)/lib/
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base
LDLIBS=$(PSBLDLIBS)

FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG).


EXEDIR=./runs


all: runsd d_matgen 

runsd:
	(if test ! -d runs ; then mkdir runs; fi)

psb_d_xyz_impl.o d_matgen.o: psb_d_xyz_mat_mod.o 

d_matgen: d_matgen.o psb_d_xyz_mat_mod.o psb_d_xyz_impl.o   
	$(FLINK) $(LINKOPT) d_matgen.o psb_d_xyz_mat_mod.o psb_d_xyz_impl.o \
	 -o d_matgen  $(PSBLAS_LIB) $(LDLIBS) 
	/bin/cp -p $(CPUPDFLAG) d_matgen $(EXEDIR)	
#	/bin/mv d_matgen $(EXEDIR)	

check: all
	cd runs && echo 5 | ./d_matgen

clean: 
	/bin/rm -f  d_matgen.o  \
	 psb_d_xyz_mat_mod.o psb_d_xyz_impl.o *$(.mod)  
verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



