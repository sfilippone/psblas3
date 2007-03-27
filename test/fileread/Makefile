include ../../Make.inc
#
# Libraries used
#
LIBDIR=../../lib/
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base

FINCLUDES=$(FMFLAG)$(LIBDIR) $(FMFLAG).

DFOBJS=getp.o  df_sample.o   
ZFOBJS=getp.o  zf_sample.o   

EXEDIR=./runs

all: df_sample zf_sample 


df_sample: $(DFOBJS)
	$(F90LINK) $(DFOBJS) -o df_sample $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv df_sample $(EXEDIR) 
zf_sample: $(ZFOBJS)
	$(F90LINK) $(ZFOBJS) -o zf_sample $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv zf_sample $(EXEDIR) 

.f90.o:
	$(MPF90) $(F90COPT) $(FINCLUDES) $(FDEFINES) -c $<

clean: 
	/bin/rm -f $(DFOBJS) $(ZFOBJS) \
        *$(.mod) $(EXEDIR)/df_sample $(EXEDIR)/zf_sample  

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)

