include Make.inc

library:
	( [ -d lib ] || mkdir lib)
	(cd src; make lib)
	@echo "====================================="
	@echo "Compilation Succesfull."
	@echo "You can now link to ./lib/libpsblas.a"

clean: 
	(cd src; make clean)

veryclean: 
	(cd src; make veryclean)
	(cd lib; /bin/rm -f *.a *$(.mod) *$(.fh))

