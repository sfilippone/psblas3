include Make.inc

all: libd based precd kryld utild cbindd
	@echo "====================================="
	@echo "PSBLAS libraries Compilation Successful."

based: libd
precd utild: based
kryld: precd based

cbindd: precd kryld utild 

libd:
	(if test ! -d lib ; then mkdir lib; fi)
	(if test ! -d include ; then mkdir include; fi; $(INSTALL_DATA) Make.inc  include/Make.inc.psblas)
	(if test ! -d modules ; then mkdir modules; fi;)	
based:
	$(MAKE) -C  base lib
precd:
	$(MAKE) -C prec lib
kryld:
	$(MAKE) -C krylov lib
utild:
	$(MAKE) -C util lib 
cbindd:
	$(MAKE) -C cbind lib 

install: all
	mkdir -p  $(INSTALL_INCLUDEDIR) &&\
	   $(INSTALL_DATA) Make.inc  $(INSTALL_INCLUDEDIR)/Make.inc.psblas
	mkdir -p  $(INSTALL_LIBDIR) &&\
	   $(INSTALL_DATA) lib/*.a  $(INSTALL_LIBDIR)
	mkdir -p  $(INSTALL_MODULESDIR) && \
	   $(INSTALL_DATA) modules/*$(.mod) $(INSTALL_MODULESDIR)
	mkdir -p  $(INSTALL_INCLUDEDIR) && \
	   $(INSTALL_DATA) include/*.h $(INSTALL_INCLUDEDIR)
	mkdir -p  $(INSTALL_DOCSDIR) && \
	   /bin/cp -fr docs/*pdf docs/html $(INSTALL_DOCSDIR)
	mkdir -p  $(INSTALL_DOCSDIR) && \
	   $(INSTALL_DATA) README.md LICENSE  $(INSTALL_DOCSDIR)
	mkdir -p  $(INSTALL_SAMPLESDIR) && \
	     /bin/cp -fr test/pargen test/fileread  $(INSTALL_SAMPLESDIR) && \
	     mkdir -p  $(INSTALL_SAMPLESDIR)/cbind && /bin/cp -fr cbind/test/pargen/* $(INSTALL_SAMPLESDIR)/cbind
clean: 
	$(MAKE) -C base clean
	$(MAKE) -C prec clean 
	$(MAKE) -C krylov clean
	$(MAKE) -C util clean
	$(MAKE) -C cbind clean

check: all
	make check -C test/serial

cleanlib:
	(cd lib; /bin/rm -f *.a *$(.mod) *$(.fh) *.h)
	(cd include; /bin/rm -f *.a *$(.mod) *$(.fh) *.h)
	(cd modules; /bin/rm -f *.a *$(.mod) *$(.fh) *.h)	

veryclean: cleanlib
	cd base && $(MAKE) veryclean
	cd prec && $(MAKE) veryclean 
	cd krylov && $(MAKE) veryclean
	cd util && $(MAKE) veryclean
	cd cbind && $(MAKE) veryclean
	cd test/fileread && $(MAKE) clean
	cd test/pargen && $(MAKE) clean
	cd test/util && $(MAKE) clean

