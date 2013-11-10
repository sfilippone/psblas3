#
# Libraries used
#
BASEDIR=../..
LIBDIR=$(BASEDIR)/lib/
INCDIR=$(BASEDIR)/include/
include $(INCDIR)/Make.inc.psblas2
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base
LDLIBS=$(PSBLDLIBS)

FINCLUDES=$(FMFLAG)$(INCDIR) $(FMFLAG).


DTOBJS=d_file_spmv.o   	
STOBJS=s_file_spmv.o   	
DPGOBJS=pdgenspmv.o
EXEDIR=./runs

all: d_file_spmv s_file_spmv pdgenspmv


d_file_spmv: $(DTOBJS)
	$(F90LINK) $(LOPT) $(DTOBJS) -o d_file_spmv $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv d_file_spmv $(EXEDIR) 

pdgenspmv: $(DPGOBJS)
	$(F90LINK) $(LOPT) $(DPGOBJS) -o pdgenspmv $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv pdgenspmv $(EXEDIR) 


s_file_spmv: $(STOBJS)
	$(F90LINK) $(LOPT) $(STOBJS) -o s_file_spmv $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv s_file_spmv $(EXEDIR) 



clean: 
	/bin/rm -f $(DBOBJSS) $(DBOBJS) $(DTOBJS) $(STOBJS)

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)

