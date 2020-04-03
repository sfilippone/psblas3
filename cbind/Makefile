include ../Make.inc

HERE=.
LIBDIR=../lib
INCDIR=../include
MODDIR=../modules/
LIBNAME=$(CBINDLIBNAME)

lib: based precd krylovd utild
	/bin/cp -p $(CPUPDFLAG) $(HERE)/$(LIBNAME) $(LIBDIR)
	/bin/cp -p $(CPUPDFLAG) *.h $(INCDIR)
	/bin/cp -p $(CPUPDFLAG) *$(.mod) $(MODDIR)


based:
	cd base && $(MAKE) lib LIBNAME=$(LIBNAME)
precd: based
	cd prec && $(MAKE) lib LIBNAME=$(LIBNAME)
krylovd: based precd
	cd krylov && $(MAKE) lib LIBNAME=$(LIBNAME)
utild: based
	cd util && $(MAKE) lib LIBNAME=$(LIBNAME)

clean:
	cd base && $(MAKE) clean
	cd prec && $(MAKE) clean
	cd krylov && $(MAKE) clean
	cd util && $(MAKE) clean

veryclean: clean
	cd test/pargen && $(MAKE) clean
	/bin/rm -f $(HERE)/$(LIBNAME) $(LIBMOD) *$(.mod) *.h
