TOP ?= ../../..
PSBLAS_DIR ?= $(TOP)
ifneq ($(wildcard $(PSBLAS_DIR)/Make.inc),)
PSBLAS_MAKE_INC = $(PSBLAS_DIR)/Make.inc
else
PSBLAS_MAKE_INC = $(PSBLAS_DIR)/include/Make.inc.psblas
endif
include $(PSBLAS_MAKE_INC)

LIBDIR = $(PSBLAS_DIR)/lib
INCLUDEDIR = $(PSBLAS_DIR)/include
PSBC_LIBS = -L$(LIBDIR) -lpsb_cbind
PSB_LIBS = -lpsb_util -lpsb_linsolve -lpsb_ext -lpsb_prec -lpsb_base
CINCLUDES = -I$(INCLUDEDIR)
PSBLAS_ARCHIVES = $(addprefix $(LIBDIR)/,libpsb_cbind.a libpsb_util.a \
  libpsb_linsolve.a libpsb_prec.a libpsb_ext.a libpsb_base.a)
PSBLAS_ARCHIVES += $(wildcard $(LIBDIR)/libpsb_cuda.a)
PSBLAS_HEADERS = $(addprefix $(INCLUDEDIR)/,psb_base_cbind.h psb_c_base.h \
  psb_c_sbase.h psb_c_dbase.h psb_c_cbase.h psb_c_zbase.h \
  psb_c_scomm.h psb_c_dcomm.h psb_c_ccomm.h psb_c_zcomm.h \
  psb_config.h psb_types.h)

MPIEXEC ?= mpiexec
MPIEXEC_NUMPROC_FLAG ?= -n
MPIEXEC_FLAGS ?=
MPIEXEC_PREFLAGS ?= $(MPIEXEC_FLAGS)
MPIEXEC_POSTFLAGS ?=
NP ?= 2
PROGRAMS = d_sparse_ops s_sparse_ops c_sparse_ops z_sparse_ops

.PHONY: all check clean
all: $(PROGRAMS)

$(PROGRAMS:=.o): Makefile $(PSBLAS_HEADERS) $(PSBLAS_MAKE_INC)

d_sparse_ops.o: sparse_ops.c
	$(MPCC) $(CCOPT) $(CDEFINES) $(CINCLUDES) -std=c99 -c $< -o $@

s_sparse_ops.o: sparse_ops.c
	$(MPCC) $(CCOPT) $(CDEFINES) $(CINCLUDES) -std=c99 -DTEST_SINGLE_REAL -c $< -o $@

c_sparse_ops.o: sparse_ops.c
	$(MPCC) $(CCOPT) $(CDEFINES) $(CINCLUDES) -std=c99 -DTEST_SINGLE_COMPLEX -c $< -o $@

z_sparse_ops.o: sparse_ops.c
	$(MPCC) $(CCOPT) $(CDEFINES) $(CINCLUDES) -std=c99 -DTEST_DOUBLE_COMPLEX -c $< -o $@

$(PROGRAMS): %: %.o $(PSBLAS_ARCHIVES)
	$(FLINK) $< -o $@ $(PSBC_LIBS) $(LPSB_CUDA) $(PSB_LIBS) $(PSBGPULDLIBS) -lm

check: all
	@set -e; for test in $(PROGRAMS); do \
	  $(MPIEXEC) $(MPIEXEC_PREFLAGS) $(MPIEXEC_NUMPROC_FLAG) $(NP) ./$$test $(MPIEXEC_POSTFLAGS); \
	done

clean:
	rm -f $(PROGRAMS) $(PROGRAMS:=.o)
