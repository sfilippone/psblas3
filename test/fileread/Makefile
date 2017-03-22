BASEDIR=../..
INCDIR=$(BASEDIR)/include/
include $(INCDIR)/Make.inc.psblas
#
# Libraries used
#
LIBDIR=$(BASEDIR)/lib/
PSBLAS_LIB= -L$(LIBDIR) -L/opencoarray-1.7.4 -lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base -lcaf_mpi
LDLIBS=$(PSBLDLIBS)

FINCLUDES=$(FMFLAG)$(INCDIR) $(FMFLAG).

SFOBJS=getp.o  sf_sample.o   
DFOBJS=getp.o  df_sample.o   
CFOBJS=getp.o  cf_sample.o   
ZFOBJS=getp.o  zf_sample.o   

EXEDIR=./runs

all: sf_sample df_sample cf_sample zf_sample 

sf_sample.o df_sample.o cf_sample.o zf_sample.o: getp.o

sf_sample: $(SFOBJS)
	$(F90LINK) $(LOPT) $(SFOBJS) -o sf_sample $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv sf_sample $(EXEDIR) 
df_sample: $(DFOBJS)
	$(F90LINK) $(LOPT) $(DFOBJS) -o df_sample $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv df_sample $(EXEDIR) 
cf_sample: $(CFOBJS)
	$(F90LINK) $(LOPT) $(CFOBJS) -o cf_sample $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv cf_sample $(EXEDIR) 
zf_sample: $(ZFOBJS)
	$(F90LINK) $(LOPT) $(ZFOBJS) -o zf_sample $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv zf_sample $(EXEDIR) 

.f90.o:
	$(MPF90) $(F90COPT) $(FINCLUDES) $(FDEFINES) -c $<

clean: 
	/bin/rm -f $(DFOBJS) $(ZFOBJS) $(SFOBJS) $(CFOBJS)\
        *$(.mod) $(EXEDIR)/df_sample $(EXEDIR)/zf_sample  

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)

