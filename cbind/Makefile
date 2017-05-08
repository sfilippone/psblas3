include ../Make.inc

HERE=.
LIBDIR=../lib
INCDIR=../include
LIBNAME=$(CBINDLIBNAME)

lib: based precd
	/bin/cp -p $(CPUPDFLAG) $(HERE)/$(LIBNAME) $(LIBDIR)
	/bin/cp -p $(CPUPDFLAG) *$(.mod) *.h $(INCDIR)



based:
	cd base && $(MAKE) lib LIBNAME=$(LIBNAME) 
precd: based
	cd prec && $(MAKE) lib LIBNAME=$(LIBNAME) 


clean: 
	cd base && $(MAKE) clean
	cd prec && $(MAKE) clean


veryclean: clean 
	/bin/rm -f $(HERE)/$(LIBNAME) $(LIBMOD) *$(.mod) *.h
