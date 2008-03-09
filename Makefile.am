include Make.inc

all: library 

library:
	(cd base; make lib)
	(cd prec; make lib )
	(cd krylov; make lib)
	(cd util; make lib )
	@echo "====================================="
	@echo "PSBLAS libraries Compilation Successful."


install: 
	($(INSTALL) -d $(INSTALL_DIR)/lib &&\
	   $(INSTALL_DATA) lib/*.a  $(INSTALL_DIR)/lib)
	($(INSTALL) -d $(INSTALL_DIR)/include && \
	   $(INSTALL_DATA) lib/*$(.mod) $(INSTALL_DIR)/include)
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

