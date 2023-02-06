include ../Make.inc

HERE=.
LIBDIR=../lib
INCDIR=../include
MODDIR=../modules/
LIBNAME=$(CBINDLIBNAME)

objs: based precd krylovd utild
	/bin/cp -p $(CPUPDFLAG) *.h $(INCDIR)
	/bin/cp -p $(CPUPDFLAG) *$(.mod) $(MODDIR)
lib: objs 
	cd base && $(MAKE) lib LIBNAME=$(LIBNAME)
	cd prec && $(MAKE) lib LIBNAME=$(LIBNAME)
	cd krylov && $(MAKE) lib LIBNAME=$(LIBNAME)
	cd util && $(MAKE) lib LIBNAME=$(LIBNAME)
	/bin/cp -p $(CPUPDFLAG) $(HERE)/$(LIBNAME) $(LIBDIR)


based:
	cd base && $(MAKE) objs LIBNAME=$(LIBNAME)
precd: based
	cd prec && $(MAKE) objs LIBNAME=$(LIBNAME)
krylovd: based precd
	cd krylov && $(MAKE) objs LIBNAME=$(LIBNAME)
utild: based
	cd util && $(MAKE) objs LIBNAME=$(LIBNAME)

clean:
	cd base && $(MAKE) clean
	cd prec && $(MAKE) clean
	cd krylov && $(MAKE) clean
	cd util && $(MAKE) clean

veryclean: clean
	cd test/pargen && $(MAKE) clean
	/bin/rm -f $(HERE)/$(LIBNAME) $(LIBMOD) *$(.mod) *.h
