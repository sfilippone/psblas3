include ../../Make.inc
#
# Libraries used
#
LIBDIR=../../lib/
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base
LDLIBS=$(PSBLDLIBS)

FINCLUDES=$(FMFLAG)$(LIBDIR) $(FMFLAG).

DBOBJS=dpsbench.o   	
DTOBJS=d_file_spmv.o   	
STOBJS=s_file_spmv.o   	
DPGOBJS=pdgenspmv.o
EXEDIR=./runs

all: d_file_spmv s_file_spmv pdgenspmv


dpsbench: $(DBOBJS)
	$(F90LINK) $(LOPT) $(DBOBJS) -o dpsbench $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv dpsbench $(EXEDIR) 

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

