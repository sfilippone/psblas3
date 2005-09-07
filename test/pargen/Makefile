include ../../Make.inc
#
# Libraries used
#
LIBDIR=../../LIB/
PSBLAS_LIB= -L$(LIBDIR) -lpsblas
SPARKER_LIB= -L$(LIBDIR) -lsparker 
ZSPARKER_LIB= -L$(LIBDIR) 
BLAS90LIB=-L$(LIBDIR) -lpsblas90
METHD90LIB=-L$(LIBDIR) 
TOOLS90LIB=-L$(LIBDIR) 

#
# Compilers and such
#
CCOPT= -g
INCDIRS=-I$(LIBDIR)


EXEDIR=./RUNS


LINKOPT=$(F90COPT)

ppde90log: ppde90log.o  part_block.o 
	$(F90LINK) $(LINKOPT) ppde90log.o part_block.o  -o ppde90log\
	$(METHD90LIB) $(TOOLS90LIB) $(BLAS90LIB) \
	$(PSBLAS_LIB) $(SPARKER_LIB) $(BLAS)\
	 $(BLACS) -llmpe -lmpe
	/bin/mv ppde90log $(EXEDIR)

ppde90: ppde90.o  part_block.o 
	$(F90LINK) $(LINKOPT) ppde90.o part_block.o  -o ppde90\
	$(METHD90LIB) $(TOOLS90LIB) $(BLAS90LIB) \
	$(PSBLAS_LIB) $(SPARKER_LIB) $(BLAS)\
	 $(BLACS) 
	/bin/mv ppde90 $(EXEDIR)
ppde90s: ppde90s.o  part_block.o 
	$(F90LINK) $(LINKOPT) ppde90s.o part_block.o  -o ppde90s\
	$(METHD90LIB) $(TOOLS90LIB) $(BLAS90LIB) \
	$(PSBLAS_LIB) $(SPARKER_LIB) $(BLAS)\
	 $(BLACS) 
	/bin/mv ppde90s $(EXEDIR)

ppde90.o: $(MODS)

pp2d: pp2d.o  part_block.o 
	$(F90LINK) $(LINKOPT) pp2d.o part_block.o  -o pp2d\
	$(METHD90LIB) $(TOOLS90LIB) $(BLAS90LIB) \
	$(METHD_LIB) $(PSBLAS_LIB) $(SPARKER_LIB) $(BLAS)\
	 $(BLACS) 
	/bin/mv pp2d $(EXEDIR)

pp2ds: pp2ds.o  part_block.o 
	$(F90LINK) $(LINKOPT) pp2ds.o part_block.o  -o pp2ds\
	$(METHD90LIB) $(TOOLS90LIB) $(BLAS90LIB) \
	$(PSBLAS_LIB) $(SPARKER_LIB) $(BLAS)\
	 $(BLACS) 
	/bin/mv pp2ds $(EXEDIR)

.f90.o:
	$(MPF90) $(F90COPT) $(INCDIRS) -c $<


clean: 
	/bin/rm -f ppde90.o  pp2d.o part_block.o $(EXEDIR)/ppde90 $(EXEDIR)/pp2d
verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make lib)



