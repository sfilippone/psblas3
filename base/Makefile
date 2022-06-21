include ../Make.inc

HERE=.
LIBDIR=../lib
INCDIR=../include
MODDIR=../modules
LIBNAME=$(BASELIBNAME)

objs:  mods sr cm in pb tl

lib:  objs
	$(MAKE) -C modules lib LIBNAME=$(BASELIBNAME) F90="$(MPF90)" F90COPT="$(F90COPT) $(MPI_OPT)"
	$(MAKE) -C serial lib LIBNAME=$(BASELIBNAME)
	$(MAKE) -C comm lib LIBNAME=$(BASELIBNAME)
	$(MAKE) -C internals lib LIBNAME=$(BASELIBNAME)
	$(MAKE) -C psblas lib LIBNAME=$(BASELIBNAME)
	$(MAKE) -C tools lib LIBNAME=$(BASELIBNAME)
	/bin/cp -p $(CPUPDFLAG) $(HERE)/$(LIBNAME) $(LIBDIR)

sr cm in pb tl: mods

mods:
	$(MAKE) -C modules objs LIBNAME=$(BASELIBNAME) F90="$(MPF90)" F90COPT="$(F90COPT) $(MPI_OPT)"
sr:
	$(MAKE) -C serial objs LIBNAME=$(BASELIBNAME)
cm:
	$(MAKE) -C comm objs LIBNAME=$(BASELIBNAME)
in:
	$(MAKE) -C internals objs LIBNAME=$(BASELIBNAME)
pb:
	$(MAKE) -C psblas objs LIBNAME=$(BASELIBNAME)
tl:
	$(MAKE) -C tools objs LIBNAME=$(BASELIBNAME)

clean: 
	($(MAKE) -C modules clean)
	($(MAKE) -C comm clean)
	($(MAKE) -C internals clean)
	($(MAKE) -C tools clean)
	($(MAKE) -C serial clean)
	($(MAKE) -C psblas clean)

veryclean: clean 
	/bin/rm -f $(HERE)/$(LIBNAME) $(LIBMOD) *$(.mod)
