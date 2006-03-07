include Make.inc

lib:
	( [ -d lib ] || mkdir lib)
	(cd src; make lib)
	@echo "====================================="
	@echo "Compilation Successful."
	@echo "You can now link to ./lib/libpsblas.a"

clean: 
	(cd src; make clean)

cleanlib:
	(cd lib; /bin/rm -f *.a *$(.mod) *$(.fh))
veryclean: cleanlib
	(cd src; make veryclean)

.PHONY: lib
