INSTALLDIR=../..
INCDIR=$(INSTALLDIR)/include
MODDIR=$(INSTALLDIR)/modules/
include $(INCDIR)/Make.inc.psblas
#
# Libraries used
LIBDIR=$(INSTALLDIR)/lib
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_linsolve -lpsb_prec -lpsb_ext -lpsb_base
LDLIBS= $(PSBLDLIBS)
#
# Compilers and such
#
CCOPT= -g
FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG).


EXEDIR=./runs

all: runsd psb_d_nest_glob_test psb_d_nest_rect_test psb_d_nest_cg_test psb_d_nest_builder_test

runsd:
	(if test ! -d runs ; then mkdir runs; fi)

psb_d_nest_glob_test: psb_d_nest_glob_test.o
	$(FLINK) psb_d_nest_glob_test.o -o psb_d_nest_glob_test $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_d_nest_glob_test $(EXEDIR)

psb_d_nest_rect_test: psb_d_nest_rect_test.o
	$(FLINK) psb_d_nest_rect_test.o -o psb_d_nest_rect_test $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_d_nest_rect_test $(EXEDIR)

psb_d_nest_cg_test: psb_d_nest_cg_test.o
	$(FLINK) psb_d_nest_cg_test.o -o psb_d_nest_cg_test $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_d_nest_cg_test $(EXEDIR)

psb_d_nest_builder_test: psb_d_nest_builder_test.o
	$(FLINK) psb_d_nest_builder_test.o -o psb_d_nest_builder_test $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_d_nest_builder_test $(EXEDIR)

clean:
	/bin/rm -f psb_d_nest_glob_test.o psb_d_nest_rect_test.o psb_d_nest_cg_test.o psb_d_nest_builder_test.o *$(.mod) \
	$(EXEDIR)/psb_d_nest_glob_test $(EXEDIR)/psb_d_nest_rect_test $(EXEDIR)/psb_d_nest_cg_test $(EXEDIR)/psb_d_nest_builder_test
verycleanlib:
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)
