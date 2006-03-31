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
	  mmio.o mat_dist.o  df_sample.o   part_blk2.o
ZFOBJS=partgraph.o part_block.o read_mat.o getp.o  \
	  mmio.o mat_dist.o  zf_sample.o   part_blk2.o
IOOBJS= mmio.o hbio.o

ZH2MOBJS=zhb2mm.o $(IOOBJS)
DH2MOBJS=dhb2mm.o $(IOOBJS)
DM2HOBJS=dmm2hb.o $(IOOBJS)
ZM2HOBJS=zmm2hb.o $(IOOBJS)
MMHBOBJS=zhb2mm.o dhb2mm.o dmm2hb.o zmm2hb.o

EXEDIR=./RUNS

all: df_sample zf_sample dhb2mm zhb2mm dmm2hb zmm2hb

read_mat.o: mmio.o

df_sample: $(DFOBJS)
	$(F90LINK) $(LINKOPT) $(DFOBJS) -o df_sample\
	$(PSBLAS_LIB) $(METIS_LIB) $(BLACS) $(SLU) $(UMF) $(BLAS)
	/bin/mv df_sample $(EXEDIR) 
zf_sample: $(ZFOBJS)
	$(F90LINK) $(LINKOPT) $(ZFOBJS) -o zf_sample\
	$(PSBLAS_LIB) $(METIS_LIB) $(BLACS) $(SLU) $(UMF) $(BLAS)
	/bin/mv zf_sample $(EXEDIR) 

$(MMHBOBJS): $(IOOBJS)
dhb2mm: $(DH2MOBJS)
	$(MPF90) -o dhb2mm   $(DH2MOBJS) $(PSBLAS_LIB) $(BLACS)
dmm2hb: $(DM2HOBJS)
	$(MPF90) -o dmm2hb   $(DM2HOBJS) $(PSBLAS_LIB) $(BLACS)
zhb2mm: $(ZH2MOBJS)
	$(MPF90) -o zhb2mm   $(ZH2MOBJS) $(PSBLAS_LIB) $(BLACS)
zmm2hb: $(ZM2HOBJS)
	$(MPF90) -o zmm2hb   $(ZM2HOBJS) $(PSBLAS_LIB) $(BLACS)

.f90.o:
	$(MPF90) $(F90COPT) $(INCDIRS) -c $<

clean: 
	/bin/rm -f $(DFOBJS) $(ZFOBJS) $(IOOBJS) $(MMHBOBJS) \
        *$(.mod) $(EXEDIR)/df_sample $(EXEDIR)/zf_sample 

lib:
	(cd ../../; make lib)
verycleanlib:
	(cd ../../; make veryclean)

