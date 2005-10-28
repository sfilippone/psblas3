include ../../Make.inc


LIBDIR=../../lib/

MPFOBJS=psb_dcslu.o psb_dbldaggrmat.o 
F90OBJS= psb_dcsrsetup.o psb_dprec.o \
	psb_dprecbld.o gps.o psb_dprecfree.o psb_dprecset.o \
	psb_dgenaggrmap.o  psb_dsplu.o $(MPFOBJS) 
#dcoocp.o dcoocpadd.o dcoofact.o dcoolu.o dcooluadd.o\

COBJS=fort_slu_impl.o fort_umf_impl.o
INCDIRS=-I. -I.. -I$(LIBDIR) 

OBJS=$(F90OBJS) $(COBJS)

lib:  mpobjs $(OBJS) 
	$(AR) $(LIBDIR)/$(LIBNAME) $(MPFOBJS) $(OBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)

mpobjs: 
	(make $(MPFOBJS) F90="$(MPF90)" F90COPT="$(F90COPT)")

veryclean: clean
	/bin/rm -f $(LIBNAME)

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

veryclean: clean
