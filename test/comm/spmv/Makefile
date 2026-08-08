INSTALLDIR=../../..
INCDIR=$(INSTALLDIR)/include/
MODDIR=$(INSTALLDIR)/modules/
include $(INCDIR)/Make.inc.psblas
#
# Libraries used
#
LIBDIR=$(INSTALLDIR)/lib/
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_linsolve -lpsb_prec -lpsb_cuda -lpsb_ext -lpsb_base
LDLIBS=$(PSBGPULDLIBS)

FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG).

TOBJS=psb_spmv_test.o

# Per-type drivers generated from the template (X_spmv_test.p90).
#
# Built CPU-only, without -DPSB_HAVE_CUDA and without -lpsb_cuda: the CUDA
# layer has no complex support upstream (psb_?_cuda_vect_mod declares an
# external cdot/zdot, which are not BLAS symbols, and cuda/Makefile never
# compiles the c/z hdiag modules), so c/z cannot be linked against it.
# Building all four the same way keeps the cross-type comparison honest.
# psb_ext stays in: it provides the CPU storage formats (ELL, HLL, DIA, HDIA).
TYPED=psb_s_spmv_test psb_d_spmv_test psb_c_spmv_test psb_z_spmv_test
TYPEDOBJS=$(addsuffix .o,$(TYPED))
CPU_PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_linsolve -lpsb_prec -lpsb_ext -lpsb_base
CPU_LDLIBS=$(PSBLDLIBS)

# GPU variants of the same drivers, s/d only for the reason above.
# Same sources, compiled with -DPSB_HAVE_CUDA and linked against libpsb_cuda,
# into differently named objects and executables so the two builds coexist.
# The source declares a module (psb_?spmv_overlap_test), so the GPU build sends
# its .mod files to GPUMODDIR: without that the two variants of the same type
# would overwrite each other's module file.
TYPED_GPU=psb_s_spmv_test_gpu psb_d_spmv_test_gpu
TYPED_GPUOBJS=$(addsuffix .o,$(TYPED_GPU))
GPUMODDIR=./gpumod

NP ?= 4

EXEDIR=./runs

all: runsd psb_spmv_kernel $(TYPED)

gpu: runsd $(TYPED_GPU)

both: all gpu

runsd:
	(if test ! -d runs ; then mkdir runs; fi)

psb_spmv_test.o: psb_spmv_test.f90
	$(FC) $(FCOPT) -cpp $(FINCLUDES) $(FDEFINES) $(FCUDEFINES) -c $< -o $@

$(TYPEDOBJS): %.o: %.f90
	$(FC) $(FCOPT) -cpp $(FINCLUDES) $(FDEFINES) -c $< -o $@

$(TYPED_GPUOBJS): %_gpu.o: %.f90
	(if test ! -d $(GPUMODDIR) ; then mkdir $(GPUMODDIR); fi)
	$(FC) $(FCOPT) -cpp -J$(GPUMODDIR) $(FINCLUDES) $(FDEFINES) $(FCUDEFINES) -c $< -o $@

psb_spmv_kernel: $(TOBJS)
	$(FLINK) $(LOPT) $(TOBJS) -o psb_spmv_kernel $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_spmv_kernel $(EXEDIR)

$(TYPED): %: %.o
	$(FLINK) $(LOPT) $< -o $@ $(CPU_PSBLAS_LIB) $(CPU_LDLIBS)
	/bin/mv $@ $(EXEDIR)

$(TYPED_GPU): %: %.o
	$(FLINK) $(LOPT) $< -o $@ $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv $@ $(EXEDIR)

run-typed: all
	for t in $(TYPED); do \
	  echo "=== $$t"; \
	  mpirun -np $(NP) $(EXEDIR)/$$t || exit 1; \
	done

clean:
	/bin/rm -f $(TOBJS) $(TYPEDOBJS) $(TYPED_GPUOBJS) $(TOBJS_API) *$(.mod) \
	  $(EXEDIR)/spmv_overlap $(EXEDIR)/psb_spmv_kernel \
	  $(addprefix $(EXEDIR)/,$(TYPED)) $(addprefix $(EXEDIR)/,$(TYPED_GPU))
	/bin/rm -rf $(GPUMODDIR)

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)