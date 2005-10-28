include ../../Make.inc
#
# Libraries used
#
LIBDIR=../../lib/
PSBLAS_LIB= -L$(LIBDIR) -lpsblas

#
# Compilers and such
#
CCOPT= -g
INCDIRS=-I$(LIBDIR)


EXEDIR=./RUNS


LINKOPT=$(F90COPT) 

ppde90: ppde90.o  part_block.o 
	$(F90LINK) $(LINKOPT) ppde90.o part_block.o  -o ppde90\
	$(PSBLAS_LIB) $(BLACS) $(SLU) $(UMF) $(BLAS)
	/bin/mv ppde90 $(EXEDIR)

ppde90.o: $(MODS)

.f90.o:
	$(MPF90) $(F90COPT) $(INCDIRS) -c $<


clean: 
	/bin/rm -f ppde90.o  pp2d.o part_block.o $(EXEDIR)/ppde90 $(EXEDIR)/pp2d
verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



