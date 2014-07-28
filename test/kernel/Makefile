#
# Libraries used
#
BASEDIR=../..
INCDIR=$(BASEDIR)/include/
include $(INCDIR)/Make.inc.psblas
LIBDIR=$(BASEDIR)/lib/
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base 
LDLIBS=$(PSBLDLIBS)

FINCLUDES=$(FMFLAG)$(INCDIR) $(FMFLAG).


DTOBJS=d_file_spmv.o   	
STOBJS=s_file_spmv.o   	
DPGOBJS=pdgenspmv.o
PWFOBJS = power_file.o
PWFROBJS = power_file_real.o
ARNOBJS =arnoldi_file.o
ARNROBJS =arnoldi_file_real.o
SIROBJS=shift_invert_real.o
EXEDIR=./runs

all: d_file_spmv s_file_spmv pdgenspmv power_file power_file_real 


d_file_spmv: $(DTOBJS)
	$(F90LINK) $(LOPT) $(DTOBJS) -o d_file_spmv $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv d_file_spmv $(EXEDIR) 

pdgenspmv: $(DPGOBJS)
	$(F90LINK) $(LOPT) $(DPGOBJS) -o pdgenspmv $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv pdgenspmv $(EXEDIR) 

s_file_spmv: $(STOBJS)
	$(F90LINK) $(LOPT) $(STOBJS) -o s_file_spmv $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv s_file_spmv $(EXEDIR) 

power_file: $(PWFOBJS)
	$(F90LINK) $(LOPT) $(PWFOBJS) -o power_file $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv power_file $(EXEDIR) 

power_file_real: $(PWFROBJS)
	$(F90LINK) $(LOPT) $(PWFROBJS) -o power_file_real $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv power_file_real $(EXEDIR) 

arnoldi_file: $(ARNOBJS)
	$(F90LINK) $(LOPT) $(ARNOBJS) -o arnoldi_file $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv arnoldi_file $(EXEDIR) 

arnoldi_file_real: $(ARNROBJS)
	$(F90LINK) $(LOPT) $(ARNROBJS) -o arnoldi_file_real $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv arnoldi_file_real $(EXEDIR) 

shift_invert_real: $(SIROBJS)
	$(F90LINK) $(LOPT) $(SIROBJS) -o shift_invert_real $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv shift_invert_real $(EXEDIR) 

analyse: analyse.o
	$(F90LINK) $(LOPT) analyse.o -o analyse $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv analyse /server/tamestoy/TEMP/MaxEigenAdj/ 

laplacian: laplacian.o
	$(F90LINK) $(LOPT) laplacian.o -o laplacian $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv laplacian $(EXEDIR) 

test_chseqr: test_chseqr.o
	$(F90LINK) $(LOPT) test_chseqr.o -o test_chseqr $(PSBLAS_LIB) $(LDLIBS)  

clean: 
	/bin/rm -f *.o
	#/bin/rm -f $(DBOBJSS) $(DBOBJS) $(PWOBJS) $(PWFOBJS) $(DTOBJS) $(STOBJS) $(ARNOBJS)

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)
