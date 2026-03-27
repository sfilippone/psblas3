include ../Make.inc


HERE=.
LIBDIR=../lib
INCDIR=../include
MODDIR=../modules

MODOBJS= psb_base_linsolve_conv_mod.o \
	psb_s_linsolve_conv_mod.o psb_c_linsolve_conv_mod.o \
	psb_d_linsolve_conv_mod.o psb_z_linsolve_conv_mod.o \
	psb_linsolve_mod.o \
	psb_pMPK_mod.o psb_eigsolve_mod.o

OBJS=$(MODOBJS)

LOCAL_MODS=$(MODOBJS:.o=$(.mod))
LIBNAME=$(METHDLIBNAME)	

FINCLUDES=$(FMFLAG). $(FMFLAG)$(MODDIR) 

objs: mods $(OBJS)  impld

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

psb_s_linsolve_conv_mod.o psb_c_linsolve_conv_mod.o psb_d_linsolve_conv_mod.o psb_z_linsolve_conv_mod.o:  psb_base_linsolve_conv_mod.o
psb_linsolve_conv_mod.o: psb_s_linsolve_conv_mod.o psb_c_linsolve_conv_mod.o psb_d_linsolve_conv_mod.o psb_z_linsolve_conv_mod.o 
$(F90OBJS): $(MODOBJS)
$(OBJS): $(MODDIR)/$(PRECMODNAME)$(.mod)  $(MODDIR)/$(BASEMODNAME)$(.mod)

veryclean: clean
	/bin/rm -f $(HERE)/$(LIBNAME)

iclean:
	cd impl && $(MAKE) clean
clean: iclean
	/bin/rm -f $(OBJS) *$(.mod)

