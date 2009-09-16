include ../Make.inc

HERE=.
LIBDIR=../lib
LIBNAME=$(BASELIBNAME)
LIBMOD=psb_base_mod$(.mod)
lib:
	(cd modules; make lib LIBNAME=$(BASELIBNAME) F90=$(MPF90) F90COPT="$(F90COPT) $(MPI_OPT)")
	(cd serial; make lib LIBNAME=$(BASELIBNAME))
#	(cd newserial; make lib LIBNAME=$(BASELIBNAME) )
	(cd comm; make lib LIBNAME=$(BASELIBNAME))
	(cd internals; make lib LIBNAME=$(BASELIBNAME))
	(cd tools; make lib LIBNAME=$(BASELIBNAME))
	(cd psblas; make lib LIBNAME=$(BASELIBNAME))
	/bin/cp -p $(HERE)/$(LIBNAME) $(LIBDIR)
	/bin/cp -p $(LIBMOD) *$(.mod) $(LIBDIR)

clean: 
	(cd modules; make clean)
	(cd comm; make clean)
	(cd internals; make clean)
	(cd tools; make clean)
	(cd serial; make clean)
	(cd psblas; make clean)
#	(cd newserial; make clean)

veryclean: clean 
	/bin/rm -f $(HERE)/$(LIBNAME) $(LIBMOD) *$(.mod)
