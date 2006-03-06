include ../../Make.inc
#
# Libraries used
#
LIBDIR=../../lib/
PSBLAS_LIB= -L$(LIBDIR) -lpsblas

#
#  We are using the public domain tool METIS from U. Minnesota. To get it
#  check URL http://www.cs.umn.edu:~karypis
#  
METIS_LIB = -L$(HOME)/NUMERICAL/metis-4.0 -lmetis 

INCDIRS=-I$(LIBDIR)

DFOBJS=partgraph.o part_block.o read_mat.o getp.o  \
	  mmio.o mat_dist.o  df_sample.o lowerc.o   part_blk2.o

EXEDIR=./RUNS

all: df_sample

read_mat.o: mmio.o

df_sample: $(DFOBJS)
	$(F90LINK) $(LINKOPT) $(DFOBJS) -o df_sample\
	$(PSBLAS_LIB) $(METIS_LIB) $(BLACS) $(SLU) $(UMF) $(BLAS)
	/bin/mv df_sample $(EXEDIR) 


aggr.o: mmio.o
aggr: aggr.o mmio.o
	$(F90LINK) $(LINKOPT) aggr.o mmio.o -o aggr \
	$(METHD90LIB) $(TOOLS90LIB) $(BLAS90LIB) $(METIS_LIB)\
	$(PSBLAS_LIB) $(SPARKER_LIB) $(BLAS)\
	 $(BLACS) 
	/bin/mv aggr $(EXEDIR) 

.f90.o:
	$(MPF90) $(F90COPT) $(INCDIRS) -c $<
#mmio.o: mmio.f90
#	$(MPF90) -O0        $(INCDIRS) -c $<


#$(DFLOBJS) $(ZFOBJS) $(DFOBJS):$(MODS)

clean: 
	/bin/rm -f $(FOBJS) ppde90.o  part_block.o $(DFOBJS) $(DFBOBJS) $(DFCOBJS) $(ZHOBJS) $(DAOBJS)\
        *$(.mod) $(EXEDIR)/df_sample $(EXEDIR)/zf_sample $(EXEDIR)/df_comm $(EXEDIR)/df_samplelog $(EXEDIR)/df_bench

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)

