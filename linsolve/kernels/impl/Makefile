include ../../../Make.inc

LIBDIR=../../../lib
INCDIR=../../../include
MODDIR=../../../modules
HERE=../..
OBJS=psb_d_pMPK_impl.o

LIBNAME=$(METHDLIBNAME)
COBJS=
FINCLUDES=$(FMFLAG).. $(FMFLAG)$(MODDIR)

objs:  $(OBJS)
lib:  objs
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)
	/bin/cp -p $(CPUPDFLAG) $(HERE)/$(LIBNAME) $(LIBDIR)

veryclean: clean

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

