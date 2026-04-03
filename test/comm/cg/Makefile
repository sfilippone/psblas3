INSTALLDIR=../../..
INCDIR=$(INSTALLDIR)/include/
MODDIR=$(INSTALLDIR)/modules/
include $(INCDIR)/Make.inc.psblas

LIBDIR=$(INSTALLDIR)/lib/
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_linsolve -lpsb_prec -lpsb_base
LDLIBS=$(PSBLDLIBS)

FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG).

NP ?= 4
IDIM ?= 40

PROGSRC=psb_comm_cg_test.F90
TOBJS=psb_comm_cg_test.o

EXEDIR=./runs
EXE=psb_comm_cg_test

all: runsd $(EXE)

runsd:
	(if test ! -d runs ; then mkdir runs; fi)

psb_comm_cg_test.o: $(PROGSRC)
	$(FC) $(FCOPT) $(FINCLUDES) $(FDEFINES) -c $(PROGSRC) -o $@

$(EXE): $(TOBJS)
	$(FLINK) $(LOPT) $(TOBJS) -o $(EXE) $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv $(EXE) $(EXEDIR)

run: all
	mpirun -np $(NP) $(EXEDIR)/$(EXE) $(IDIM)

clean:
	/bin/rm -f $(TOBJS) *$(.mod) $(EXEDIR)/$(EXE)
