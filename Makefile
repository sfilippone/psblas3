include Make.inc

all: dirs mods objs libd
	@echo "====================================="
	@echo "PSBLAS libraries Compilation Successful."

dirs:
	(if test ! -d lib ; then mkdir lib; fi)
	(if test ! -d include ; then mkdir include; fi; $(INSTALL_DATA) Make.inc  include/Make.inc.psblas)
	(if test ! -d modules ; then mkdir modules; fi;)	

mods: basemods utilmods precmods linslvmods cbindmods extmods $(CUDAMODS) $(OACCMODS)

basemods:
	$(MAKE) -C base mods
precmods: basemods
	$(MAKE) -C prec mods
linslvmods: precmods
	$(MAKE) -C linsolve mods
utilmods: basemods
	$(MAKE) -C util mods 
cbindmods: basemods precmods linslvmods utilmods extmods $(CUDAMODS)
	$(MAKE) -C cbind objs
extmods: basemods  
	$(MAKE) -C ext mods
cudamods: extmods  
	$(MAKE) -C cuda mods
oaccmods: extmods  
	$(MAKE) -C openacc mods 


objs: mods based precd linslvd utild cbindd extd  $(CUDAD) $(OACCD)
based: basemods
precd: precmods
utild: utilmods	
linslvd: linslvmods
extd:  extmods
cudad:  cudamods
oaccd:  oaccmods
cbindd: cbindmods

libd: based precd linslvd utild cbindd extd $(CUDALD) $(OACCLD)
	$(MAKE) -C base lib
	$(MAKE) -C prec lib
	$(MAKE) -C linsolve lib
	$(MAKE) -C util lib 
	$(MAKE) -C cbind lib
	$(MAKE) -C ext lib
cudald:  cudad
	$(MAKE) -C cuda lib
oaccld:  oaccd
	$(MAKE) -C openacc lib


based:
	$(MAKE) -C base objs
precd:
	$(MAKE) -C prec objs
linslvd:
	$(MAKE) -C linsolve objs
utild:
	$(MAKE) -C util objs 
cbindd:
	$(MAKE) -C cbind objs 
extd:   
	$(MAKE) -C ext objs
cudad:   cudamods
	$(MAKE) -C cuda objs
oaccd:   oaccmods
	$(MAKE) -C openacc objs 


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
	     /bin/cp -fr test/pdegen test/fileread  $(INSTALL_SAMPLESDIR) && \
	     mkdir -p  $(INSTALL_SAMPLESDIR)/cbind && /bin/cp -fr cbind/test/pdegen/* $(INSTALL_SAMPLESDIR)/cbind
clean: cleanlib
	$(MAKE) -C base veryclean
	$(MAKE) -C prec veryclean 
	$(MAKE) -C linsolve veryclean
	$(MAKE) -C util veryclean
	$(MAKE) -C cbind veryclean
	$(MAKE) -C ext veryclean
	$(MAKE) -C cuda veryclean
	$(MAKE) -C openacc veryclean
cleantest:
	cd test/fileread && $(MAKE) clean
	cd test/pdegen && $(MAKE) clean
	cd test/util && $(MAKE) clean

cleanlib:
	(cd lib; /bin/rm -f *.a *$(.mod) *$(.fh) *.h)
	(cd include; /bin/rm -f *.a *$(.mod) *$(.fh) *.h)
	(cd modules; /bin/rm -f *.a *$(.mod) *$(.fh) *.h)

distclean: clean
	/bin/rm -f Make.inc  util/psb_metis_int.h base/modules/psb_config.h 

check: all
	make check -C test/serial


