INSTALLDIR=../../..
INCDIR=$(INSTALLDIR)/include/
MODDIR=$(INSTALLDIR)/modules/
include $(INCDIR)/Make.inc.psblas

LIBDIR=$(INSTALLDIR)/lib/
PSBLAS_LIB= -L$(LIBDIR) $(LCUDA) -lpsb_util -lpsb_linsolve -lpsb_prec -lpsb_base -lpsb_cuda -lpsb_ext
LDLIBS=$(PSBGPULDLIBS) 

FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG).

NP ?= 4
IDIM ?= 40

PROGSRC=psb_comm_cg_test.F90
TOBJS=psb_comm_cg_test.o

# Per-type drivers generated from the template (X_comm_cg_test.P90).
#
# Built CPU-only, without -DPSB_HAVE_CUDA and without -lpsb_cuda: the CUDA
# layer has no complex support upstream (psb_?_cuda_vect_mod declares an
# external cdot/zdot, which are not BLAS symbols, and cuda/Makefile never
# compiles the c/z hdiag modules), so c/z cannot be linked against it.
# Building all four the same way keeps the cross-type comparison honest.
# psb_ext stays in: it provides the CPU storage formats (ELL, HLL, DIA, HDIA).
TYPED=psb_s_comm_cg_test psb_d_comm_cg_test psb_c_comm_cg_test psb_z_comm_cg_test
TYPEDOBJS=$(addsuffix .o,$(TYPED))
CPU_PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_linsolve -lpsb_prec -lpsb_ext -lpsb_base
CPU_LDLIBS=$(PSBLDLIBS)

# GPU variants of the same drivers, s/d only for the reason above. Same
# sources, compiled with -DPSB_HAVE_CUDA and linked against libpsb_cuda, into
# differently named objects and executables so the two builds coexist. These
# drivers declare no module, so no module-file separation is needed.
TYPED_GPU=psb_s_comm_cg_test_gpu psb_d_comm_cg_test_gpu
TYPED_GPUOBJS=$(addsuffix .o,$(TYPED_GPU))

EXEDIR=./runs
EXE=psb_comm_cg_test

all: runsd $(EXE) $(TYPED)

gpu: runsd $(TYPED_GPU)

both: all gpu

runsd:
	(if test ! -d runs ; then mkdir runs; fi)

psb_comm_cg_test.o: $(PROGSRC)
	$(FC) $(FCOPT) $(FINCLUDES) $(FDEFINES) $(FCUDEFINES) -c $(PROGSRC) -o $@

$(TYPEDOBJS): %.o: %.F90
	$(FC) $(FCOPT) $(FINCLUDES) $(FDEFINES) -c $< -o $@

$(TYPED_GPUOBJS): %_gpu.o: %.F90
	$(FC) $(FCOPT) $(FINCLUDES) $(FDEFINES) $(FCUDEFINES) -c $< -o $@

$(EXE): $(TOBJS)
	$(FLINK) $(LOPT) $(TOBJS) -o $(EXE) $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv $(EXE) $(EXEDIR)

$(TYPED): %: %.o
	$(FLINK) $(LOPT) $< -o $@ $(CPU_PSBLAS_LIB) $(CPU_LDLIBS)
	/bin/mv $@ $(EXEDIR)

$(TYPED_GPU): %: %.o
	$(FLINK) $(LOPT) $< -o $@ $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv $@ $(EXEDIR)

run: all
	mpirun -np $(NP) $(EXEDIR)/$(EXE) $(IDIM)

run-typed: all
	for t in $(TYPED); do \
	  echo "=== $$t"; \
	  mpirun -np $(NP) $(EXEDIR)/$$t $(IDIM) || exit 1; \
	done

clean:
	/bin/rm -f $(TOBJS) $(TYPEDOBJS) $(TYPED_GPUOBJS) *$(.mod) \
	  $(EXEDIR)/$(EXE) $(addprefix $(EXEDIR)/,$(TYPED)) \
	  $(addprefix $(EXEDIR)/,$(TYPED_GPU))
