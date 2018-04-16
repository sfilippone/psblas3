include Make.inc

all: libd based precd kryld utild 
	@echo "====================================="
	@echo "PSBLAS libraries Compilation Successful."
#cbindd
based: libd
precd utild: based
kryld: precd based

cbindd: precd kryld utild 

libd:
	(if test ! -d lib ; then mkdir lib; fi)
	(if test ! -d include ; then mkdir include; fi; $(INSTALL_DATA) Make.inc  include/Make.inc.psblas)
	(if test ! -d modules ; then mkdir modules; fi;)	
based:
	cd base && $(MAKE) lib
precd:
	cd prec && $(MAKE) lib
kryld:
	cd krylov && $(MAKE) lib
utild:
	cd util&& $(MAKE) lib 
cbindd:
	cd cbind&& $(MAKE) lib 

install: all
	$(SHELL) ./mkdir.sh  $(INSTALL_INCLUDEDIR) &&\
	   $(INSTALL_DATA) Make.inc  $(INSTALL_INCLUDEDIR)/Make.inc.psblas
	$(SHELL) ./mkdir.sh  $(INSTALL_LIBDIR) &&\
	   $(INSTALL_DATA) lib/*.a  $(INSTALL_LIBDIR)
	$(SHELL) ./mkdir.sh  $(INSTALL_MODULESDIR) && \
	   $(INSTALL_DATA) modules/*$(.mod) $(INSTALL_MODULESDIR)
	$(SHELL) ./mkdir.sh  $(INSTALL_INCLUDEDIR) && \
	   $(INSTALL_DATA) include/*.h $(INSTALL_INCLUDEDIR)
	$(SHELL) ./mkdir.sh  $(INSTALL_DOCSDIR) && \
	   /bin/cp -fr docs/*pdf docs/html $(INSTALL_DOCSDIR)
	$(SHELL) ./mkdir.sh  $(INSTALL_DOCSDIR) && \
	   $(INSTALL_DATA) README.md LICENSE  $(INSTALL_DOCSDIR)
	$(SHELL) ./mkdir.sh  $(INSTALL_SAMPLESDIR) && \
	     /bin/cp -fr test/pargen test/fileread test/kernel $(INSTALL_SAMPLESDIR) && \
	     ./mkdir.sh $(INSTALL_SAMPLESDIR)/cbind && /bin/cp -fr cbind/test/pargen/* $(INSTALL_SAMPLESDIR)/cbind
clean: 
	cd base && $(MAKE) clean
	cd prec && $(MAKE) clean 
	cd krylov && $(MAKE) clean
	cd util && $(MAKE) clean
	cd cbind && $(MAKE) clean

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

