INSTALLDIR=../../..
INCDIR=$(INSTALLDIR)/include/
MODDIR=$(INSTALLDIR)/modules/
include $(INCDIR)/Make.inc.psblas
#
# Libraries used
#
LIBDIR=$(INSTALLDIR)/lib/
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_linsolve -lpsb_prec -lpsb_base
LDLIBS=$(PSBLDLIBS)

FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG).

TOBJS=psb_comm_test.o

# Per-type drivers generated from the template (X_comm_test.P90).
TYPED=psb_s_comm_test psb_d_comm_test psb_c_comm_test psb_z_comm_test
TYPEDOBJS=$(addsuffix .o,$(TYPED))

NP ?= 4

EXEDIR=./runs

all: runsd psb_comm_test $(TYPED)

runsd:
	(if test ! -d runs ; then mkdir runs; fi)

psb_comm_test: $(TOBJS)
	$(FLINK) $(LOPT) $(TOBJS) -o psb_comm_test $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_comm_test $(EXEDIR)

$(TYPED): %: %.o
	$(FLINK) $(LOPT) $< -o $@ $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv $@ $(EXEDIR)

run: all
	for t in $(TYPED); do \
	  echo "=== $$t"; \
	  mpirun -np $(NP) $(EXEDIR)/$$t || exit 1; \
	done

clean:
	/bin/rm -f $(TOBJS) $(TYPEDOBJS) $(TOBJS_API) *$(.mod) \
	  $(EXEDIR)/psb_comm_test $(addprefix $(EXEDIR)/,$(TYPED))

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)