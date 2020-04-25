include ../Make.inc

HERE=.
LIBDIR=../lib
INCDIR=../include
MODDIR=../modules
LIBNAME=$(BASELIBNAME)

lib: mods sr cm in pb tl
	/bin/cp -p $(CPUPDFLAG) $(HERE)/$(LIBNAME) $(LIBDIR)
	/bin/cp -p $(CPUPDFLAG) *$(.mod) $(MODDIR)

sr cm in pb tl: mods

mods:
	$(MAKE) -C modules lib LIBNAME=$(BASELIBNAME) F90="$(MPF90)" F90COPT="$(F90COPT) $(MPI_OPT)"
sr:
	$(MAKE) -C serial lib LIBNAME=$(BASELIBNAME)
cm:
	$(MAKE) -C comm lib LIBNAME=$(BASELIBNAME)
in:
	$(MAKE) -C internals lib LIBNAME=$(BASELIBNAME)
pb:
	$(MAKE) -C psblas lib LIBNAME=$(BASELIBNAME)
tl:
	$(MAKE) -C tools lib LIBNAME=$(BASELIBNAME)

clean: 
	($(MAKE) -C modules clean)
	($(MAKE) -C comm clean)
	($(MAKE) -C internals clean)
	($(MAKE) -C tools clean)
	($(MAKE) -C serial clean)
	($(MAKE) -C psblas clean)

veryclean: clean 
	/bin/rm -f $(HERE)/$(LIBNAME) $(LIBMOD) *$(.mod)
