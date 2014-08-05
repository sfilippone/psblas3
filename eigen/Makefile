include ../Make.inc


HERE=.
LIBDIR=../lib
INCDIR=../include

MODOBJS=psb_eigen_mod.o 
F90OBJS=psb_d_power_vect.o
OBJS=$(F90OBJS) $(MODOBJS)

LOCAL_MODS=$(MODOBJS:.o=$(.mod))

LIBNAME=$(EIGENLIBNAME)	

FINCLUDES=$(FMFLAG). $(FMFLAG)$(INCDIR) 

lib: $(OBJS) 
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)
	/bin/cp -p $(CPUPDFLAG) $(HERE)/$(LIBNAME) $(LIBDIR)
	/bin/cp -p $(CPUPDFLAG) *$(.mod) $(INCDIR)


$(F90OBJS): $(MODOBJS)
$(OBJS): $(INCDIR)/$(PRECMODNAME)$(.mod)  $(INCDIR)/$(BASEMODNAME)$(.mod)


veryclean: clean
	/bin/rm -f $(HERE)/$(LIBNAME)

clean:
	/bin/rm -f $(OBJS) *$(.mod)

