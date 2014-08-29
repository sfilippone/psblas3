#
# Libraries used
#
BASEDIR=../..
INCDIR=$(BASEDIR)/include/
include $(INCDIR)/Make.inc.psblas
LIBDIR=$(BASEDIR)/lib/
PSBLAS_LIB= -L$(LIBDIR) -lpsb_eigen -lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base 

LDLIBS=$(PSBLDLIBS)

FINCLUDES=$(FMFLAG)$(INCDIR) $(FMFLAG).

PWFROBJS = power_file_real.o 
EXEDIR=./runs

all: power_file_real maxA_and_lapl_extremums test_enea_dagostino test_realnet


power_file_real: $(PWFROBJS)
	$(F90LINK) $(LOPT) $(PWFROBJS) -o power_file_real $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv power_file_real $(EXEDIR) 

maxA_and_lapl_extremums: maxA_and_lapl_extremums.o
	$(F90LINK) $(LOPT) maxA_and_lapl_extremums.o -o maxA_and_lapl_extremums $(PSBLAS_LIB) $(LDLIBS) 
	/bin/cp maxA_and_lapl_extremums test_dagostino
	/bin/mv maxA_and_lapl_extremums $(EXEDIR) 

test_enea_dagostino:test_enea_dagostino.o
	$(F90LINK) $(LOPT) test_enea_dagostino.o -o test_enea_dagostino $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv test_enea_dagostino test_dagostino

adj_to_mm : adj_to_mm.o
	$(F90LINK) $(LOPT) adj_to_mm.o -o adj_to_mm $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv adj_to_mm $(EXEDIR)

test_realnet:test_realnet.o
	$(F90LINK) $(LOPT) test_realnet.o -o test_realnet $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv test_realnet test_dagostino

analyse_results: analyse_results.o
	$(F90LINK) $(LOPT) analyse_results.o -o analyse_results $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv analyse_results /server/tamestoy/TEMP/SCALING-BA 

clean: 
	/bin/rm -f *.o
	#/bin/rm -f $(DBOBJSS) $(DBOBJS) $(PWOBJS) $(PWFOBJS) $(DTOBJS) $(STOBJS) $(ARNOBJS)

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)
