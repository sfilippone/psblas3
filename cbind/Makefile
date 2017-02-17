include ../Make.inc

HERE=.
LIBDIR=../lib
INCDIR=../include
LIBNAME=$(CBINDLIBNAME)

lib: based 
	/bin/cp -p $(CPUPDFLAG) $(HERE)/$(LIBNAME) $(LIBDIR)
	/bin/cp -p $(CPUPDFLAG) *$(.mod) *.h $(INCDIR)



based:
	cd base && $(MAKE) lib LIBNAME=$(LIBNAME) 


clean: 
	cd base && $(MAKE) clean


veryclean: clean 
	/bin/rm -f $(HERE)/$(LIBNAME) $(LIBMOD) *$(.mod) *.h
