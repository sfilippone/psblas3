INSTALLDIR=../..
INCDIR=$(INSTALLDIR)/include/
MODDIR=$(INSTALLDIR)/modules/
include $(INCDIR)/Make.inc.psblas
#
# Libraries used
#
LIBDIR=$(INSTALLDIR)/lib/
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base
LDLIBS=$(PSBLDLIBS)

FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG).

SFOBJS=getp.o  psb_sf_sample.o   
DFOBJS=getp.o  psb_df_sample.o   
CFOBJS=getp.o  psb_cf_sample.o   
ZFOBJS=getp.o  psb_zf_sample.o   

EXEDIR=./runs

all: psb_sf_sample psb_df_sample psb_cf_sample psb_zf_sample 

psb_sf_sample.o psb_df_sample.o psb_cf_sample.o psb_zf_sample.o: getp.o

psb_sf_sample: $(SFOBJS)
	$(FLINK) $(LOPT) $(SFOBJS) -o psb_sf_sample $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv psb_sf_sample $(EXEDIR) 
psb_df_sample: $(DFOBJS)
	$(FLINK) $(LOPT) $(DFOBJS) -o psb_df_sample $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv psb_df_sample $(EXEDIR) 
psb_cf_sample: $(CFOBJS)
	$(FLINK) $(LOPT) $(CFOBJS) -o psb_cf_sample $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_cf_sample $(EXEDIR) 
psb_zf_sample: $(ZFOBJS)
	$(FLINK) $(LOPT) $(ZFOBJS) -o psb_zf_sample $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_zf_sample $(EXEDIR) 

.f90.o:
	$(MPFC) $(FCOPT) $(FINCLUDES) $(FDEFINES) -c $<

clean: 
	/bin/rm -f $(DFOBJS) $(ZFOBJS) $(SFOBJS) $(CFOBJS)\
        *$(.mod) $(EXEDIR)/psb_*f_sample 

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)

