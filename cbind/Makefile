include ../Make.inc

HERE=.
LIBDIR=../lib
INCDIR=../include
MODDIR=../modules/
LIBNAME=$(CBINDLIBNAME)
MPI_TEST_NPROCS ?= 1 2 4

.PHONY: objs lib based precd linsolved utild check clean veryclean

objs: based precd linsolved utild
	/bin/cp -p $(CPUPDFLAG) *.h $(INCDIR)
	/bin/cp -p $(CPUPDFLAG) *$(.mod) $(MODDIR)
lib: objs 
	cd base && $(MAKE) lib LIBNAME=$(LIBNAME)
	cd prec && $(MAKE) lib LIBNAME=$(LIBNAME)
	cd linsolve && $(MAKE) lib LIBNAME=$(LIBNAME)
	cd util && $(MAKE) lib LIBNAME=$(LIBNAME)
	/bin/cp -p $(CPUPDFLAG) $(HERE)/$(LIBNAME) $(LIBDIR)


based:
	cd base && $(MAKE) objs LIBNAME=$(LIBNAME)
precd: based
	cd prec && $(MAKE) objs LIBNAME=$(LIBNAME)
linsolved: based precd
	cd linsolve && $(MAKE) objs LIBNAME=$(LIBNAME)
utild: based
	cd util && $(MAKE) objs LIBNAME=$(LIBNAME)

check: lib
	@set -e; for np in $(MPI_TEST_NPROCS); do \
	  $(MAKE) -C test/sparse_ops check NP=$$np; \
	done

clean:
	cd base && $(MAKE) clean
	cd prec && $(MAKE) clean
	cd linsolve && $(MAKE) clean
	cd util && $(MAKE) clean
	$(MAKE) -C test/sparse_ops clean

veryclean: clean
	cd test/pdegen && $(MAKE) clean
	/bin/rm -f $(HERE)/$(LIBNAME) $(LIBMOD) *$(.mod) *.h
