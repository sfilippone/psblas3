include ../../Make.inc
#
# Libraries used
#
LIBDIR=../../LIB/
PSBLAS_LIB= -L$(LIBDIR) -lpsblas
SPARKER_LIB= -L$(LIBDIR) -lsparker 
BLAS90LIB=-L$(LIBDIR) -lpsblas90 $(SLU) 
#METHD90LIB=-L$(LIBDIR) -lmethd90
#TOOLS90LIB=-L$(LIBDIR) -ltools90
#PREC90LIB=-L$(LIBDIR) -lprec90



#
#  We are using the public domain tool METIS from U. Minnesota. To get it
#  check URL http://www.cs.umn.edu:~karypis
#  
METIS_LIB = -L$(HOME)/NUMERICAL/metis-4.0 -lmetis 


#
CCOPT= -g
INCDIRS=-I$(LIBDIR)


TMMOBJS=partgraph.o part_block.o read_mat.o getp.o  \
	  mmio.o mat_dist.o  testmm.o  
DFOBJS=partgraph.o part_block.o read_mat.o getp.o  \
	  mmio.o mat_dist.o  df_sample.o lowerc.o   part_blk2.o
DAOBJS=partgraph.o part_block.o  getp.o  \
	  mmio.o mat_dist.o read_mat.o d_aggr.o   part_blk2.o lowerc.o 
DFLOBJS=partgraph.o part_block.o read_mat.o getp.o \
	  mmio.o mat_dist.o df_samplelog.o  

ZFOBJS=partgraph.o part_block.o read_mat.o getp.o  \
	  mmio.o mat_dist.o zf_sample.o
DFMOBJS=partgraph.o part_block.o mmio.o read_mat.o \
	   mat_dist.o  df_samplem.o   part_blk2.o lowerc.o
DFBOBJS=partgraph.o part_block.o mmio.o read_mat.o \
	   mat_dist.o comm_info.o df_bench.o   part_blk2.o lowerc.o
DFCOBJS=partgraph.o part_block.o mmio.o read_mat.o \
	   mat_dist.o comm_info.o df_comm.o   part_blk2.o lowerc.o
        


EXEDIR=./RUNS

all: df_sample zf_sample df_samplelog testmm df_samplem

read_mat.o: mmio.o

df_sample: $(DFOBJS)
	$(F90LINK) $(LINKOPT) $(DFOBJS) -o df_sample\
	$(METHD90LIB) $(TOOLS90LIB) $(BLAS90LIB) $(METIS_LIB)\
	$(SPARKER_LIB) $(PREC90LIB) $(BLAS90LIB) $(TOOLS90LIB) $(PSBLAS_LIB) \
	$(BLAS) $(SPARKER_LIB) $(BLACS) $(BLAS)
	/bin/mv df_sample $(EXEDIR) 

d_aggr: $(DAOBJS)
	$(F90LINK) $(LINKOPT) $(DAOBJS) -o d_aggr\
	$(METHD90LIB) $(TOOLS90LIB) $(BLAS90LIB) $(METIS_LIB)\
	$(SPARKER_LIB) $(PREC90LIB) $(BLAS90LIB) $(TOOLS90LIB) $(PSBLAS_LIB) \
	$(BLAS) $(SPARKER_LIB) $(BLACS) $(BLAS)
	/bin/mv d_aggr $(EXEDIR) 


df_samplem: $(DFMOBJS)
	$(F90LINK) $(LINKOPT) $(DFMOBJS) -o df_samplem\
	$(METHD90LIB) $(TOOLS90LIB) $(BLAS90LIB) $(METIS_LIB)\
	$(SPARKER_LIB) $(PREC90LIB) $(BLAS90LIB) $(TOOLS90LIB) $(PSBLAS_LIB) \
	$(BLAS) $(SPARKER_LIB) $(BLACS) $(BLAS)
	/bin/mv df_samplem $(EXEDIR) 
df_bench: $(DFBOBJS)
	$(F90LINK) $(LINKOPT) $(DFBOBJS) -o df_bench\
	$(METHD90LIB) $(TOOLS90LIB) $(BLAS90LIB) $(METIS_LIB)\
	$(SPARKER_LIB) $(PREC90LIB) $(BLAS90LIB) $(TOOLS90LIB) $(PSBLAS_LIB) \
	$(BLAS) $(SPARKER_LIB) $(BLACS) $(BLAS)
	/bin/mv df_bench $(EXEDIR) 

df_comm: $(DFCOBJS)
	$(F90LINK) $(LINKOPT) $(DFCOBJS) -o df_comm\
	$(METHD90LIB) $(TOOLS90LIB) $(BLAS90LIB) $(METIS_LIB)\
	$(SPARKER_LIB) $(PREC90LIB) $(BLAS90LIB) $(TOOLS90LIB) $(PSBLAS_LIB) \
	$(BLAS) $(SPARKER_LIB) $(BLACS) $(BLAS)
	/bin/mv df_comm $(EXEDIR) 

testmm: $(TMMOBJS)
	$(F90LINK) $(LINKOPT) $(TMMOBJS) -o testmm\
	$(METHD90LIB) $(TOOLS90LIB) $(BLAS90LIB) $(METIS_LIB)\
	$(SPARKER_LIB) $(PREC90LIB) $(BLAS90LIB) $(TOOLS90LIB) $(PSBLAS_LIB) \
	$(BLAS) $(SPARKER_LIB) $(BLACS) $(BLAS)
	/bin/mv testmm $(EXEDIR) 

df_samplelog: $(DFLOBJS)
	$(F90LINK) $(LINKOPT) $(DFLOBJS) -o df_samplelog\
	$(METHD90LIB) $(TOOLS90LIB) $(BLAS90LIB) $(METIS_LIB)\
	$(SPARKER_LIB) $(PREC90LIB) $(BLAS90LIB) $(PSBLAS_LIB) \
	$(BLAS) $(SPARKER_LIB) $(BLACS) $(BLAS) -llmpe -lmpe
	/bin/mv df_samplelog $(EXEDIR) 

zf_sample: $(ZFOBJS)
	$(F90LINK) $(LINKOPT) $(ZFOBJS) -o zf_sample\
	$(METHD90LIB) $(TOOLS90LIB) $(BLAS90LIB) $(METIS_LIB)\
	$(PSBLAS_LIB) $(SPARKER_LIB) $(BLAS)\
	 $(BLACS) 
	/bin/mv zf_sample $(EXEDIR) 

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
	(cd ../../; make lib)
verycleanlib:
	(cd ../../; make veryclean)

