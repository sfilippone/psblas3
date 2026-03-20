include ../../Make.inc


HERE=..
LIBDIR=../../lib
INCDIR=../../include
MODDIR=../../modules

MODOBJS=psb_d_pMPK_mod.o psb_pMPK_mod.o

OBJS=$(MODOBJS)

LOCAL_MODS=$(MODOBJS:.o=$(.mod))
LIBNAME=$(METHDLIBNAME)	

FINCLUDES=$(FMFLAG). $(FMFLAG)$(MODDIR) 

objs: mods $(OBJS) impld

mods: $(MODOBJS)
	$(INSTALL) -p $(CPUPDFLAG) *$(.mod) $(MODDIR)

lib: objs
	$(MAKE) -C impl lib
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)
	/bin/cp -p $(CPUPDFLAG) $(HERE)/$(LIBNAME) $(LIBDIR)


$(OBJS): $(MODDIR)/$(BASEMODNAME)$(.mod)

impld: $(OBJS)
	$(MAKE) -C impl objs

psb_pMPK_mod.o: psb_d_pMPK_mod.o 
$(F90OBJS): $(MODOBJS)
$(OBJS): $(MODDIR)/$(PRECMODNAME)$(.mod)  $(MODDIR)/$(BASEMODNAME)$(.mod)

veryclean: clean
	/bin/rm -f $(HERE)/$(LIBNAME)

iclean:
	cd impl && $(MAKE) clean
clean: iclean
	/bin/rm -f $(OBJS) *$(.mod)

