include Make.inc

library:
	( [ -d lib ] || mkdir lib)
	(cd src; make lib)

clean: 
	(cd src; make clean)

veryclean: 
	(cd src; make veryclean)
	(cd lib; /bin/rm -f *.a *$(.mod) V*.inc *.pc *.pcl)

