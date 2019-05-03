INSTALLDIR=../..
INCDIR=$(INSTALLDIR)/include/
MODDIR=$(INSTALLDIR)/modules/
include $(INCDIR)/Make.inc.psblas
P=8
debug=
FFLAGS += -cpp -DITERATIONS=$(ITERATIONS) -DSWAPCHOICE=$(SWAPCHOICE)
FFLAGS += $(debug)
FCOPT += $(FFLAGS)
scorep_vars=SCOREP_WRAPPER_COMPILER_FLAGS="-O2" \
     SCOREP_MPI_ENABLE_GROUP=all SCOREP_TIMER=gettimeofday
scorep_dir=test-results/scorep-2019-test-$(ITERATIONS)
scorep_opt=SCOREP_EXPERIMENT_DIRECTORY=$(scorep_dir) $(scorep_vars)
#
# Libraries used
#
MPFC=scorep-mpif90
# MPFC=mpif90
LIBDIR=$(INSTALLDIR)/lib/
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base
LDLIBS=$(PSBLDLIBS)
FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG).
HALOOBJS=getp.o
EXEDIR=./runs

all: psb_halo_df_test

psb_halo_df_test: getp.o psb_halo_df_test.o
	$(scorep_opt) $(FLINK) $(LOPT) $^ -o halo_$(SWAPCHOICE)$(ITERATIONS)_iter $(PSBLAS_LIB) $(LDLIBS)

halo_test:  getp.o psb_halo_df_test.o
	$(scorep_opt) $(FLINK) $(LOPT)  $^ -o halo_test $(PSBLAS_LIB) $(LDLIBS)

clean:
	rm -f *~ *.o psb_halo_df_test1 psb_halo_df_test halo_test \
	rm -rf test-results/*
	/bin/rm -f $(DFOBJS) $(ZFOBJS) $(SFOBJS) $(CFOBJS)\
        *$(.mod) $(EXEDIR)/psb_*f_sample
cleanall: clean
	rm -f halo_output.txt
	rm -rf scorep-2019*
	rm -rf scorep-failed-*
score:
	scorep-score -r $(scorep_dir)/profile.cubex | grep -e MPI -e 'type max_buf'

run: psb_halo_df_test
	$(scorep_opt) mpirun -np $(P) ./halo_$(SWAPCHOICE)_$(ITERATIONS)_iter matrices/ck104/ck104.inp

r: halo_test
	$(scorep_opt) mpirun -np $(P) ./halo_test matrices/ck104/ck104.inp

subrun:
	$(scorep_opt) mpirun -machinefile hostfile -ppn 1 -np $(cpus) ./psb_halo_df_test_$(ITERATIONS) matrices/$(MATRIX)/$(MATRIX).inp

debug: getp.o psb_halo_df_test.o
	$(FLINK) $(debug) $(LOPT) $^ -o psb_halo_df_test $(PSBLAS_LIB) $(LDLIBS)
	mpirun -np $(P) xterm -hold -e gdb -ex run --args ./psb_halo_df_test matrices/ck104/ck104.inp

test: getp.o psb_halo_df_test1.o
	$(FLINK) $(debug) $(LOPT) $^ -o psb_halo_df_test1 $(PSBLAS_LIB) $(LDLIBS)
	mpirun -np $(P) ./psb_halo_df_test1 matrices/pde225/pde225.inp
lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)
