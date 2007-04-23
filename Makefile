include Make.inc
#PREC=../mld2p4
PREC=prec


library:
	( [ -d lib ] || mkdir lib)
	(cd base; make lib)
	(cd $(PREC); make lib )
	(cd krylov; make lib)
	(cd util; make lib )
	@echo "====================================="
	@echo "PSBLAS libraries Compilation Successful."

clean: 
	(cd base; make clean)
	(cd $(PREC); make clean )
	(cd krylov; make clean)
	(cd util; make clean)

cleanlib:
	(cd lib; /bin/rm -f *.a *$(.mod) *$(.fh))
veryclean: cleanlib
	(cd base; make veryclean)
	(cd $(PREC); make veryclean )
	(cd krylov; make veryclean)
	(cd util; make veryclean)

.PHONY: lib
