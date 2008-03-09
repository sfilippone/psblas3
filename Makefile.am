include Make.inc

all: library 

library:
	(if test ! -d lib ; then mkdir lib; fi)
	(cd base; make lib)
	(cd prec; make lib )
	(cd krylov; make lib)
	(cd util; make lib )
	@echo "====================================="
	@echo "PSBLAS libraries Compilation Successful."

install:
	($(INSTALL) -d $(INSTALL_LIBDIR) &&\
	   $(INSTALL_DATA) lib/*.a  $(INSTALL_LIBDIR))
	($(INSTALL) -d $(INSTALL_INCLUDEDIR) && \
	   $(INSTALL_DATA) lib/*$(.mod) $(INSTALL_INCLUDEDIR))
clean: 
	(cd base; make clean)
	(cd prec; make clean )
	(cd krylov; make clean)
	(cd util; make clean)

cleanlib:
	(cd lib; /bin/rm -f *.a *$(.mod) *$(.fh))
veryclean: cleanlib
	(cd base; make veryclean)
	(cd prec; make veryclean )
	(cd krylov; make veryclean)
	(cd util; make veryclean)

