BASEDIR=../..
INCDIR=$(BASEDIR)/include/
include $(INCDIR)/Make.inc.psblas
#
# Libraries used
#
LIBDIR=$(BASEDIR)/lib/
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base
LDLIBS=$(PSBLDLIBS)

FINCLUDES=$(FMFLAG)$(INCDIR) $(FMFLAG).

ZH2MOBJS=zhb2mm.o 
DH2MOBJS=dhb2mm.o 
DM2HOBJS=dmm2hb.o 
ZM2HOBJS=zmm2hb.o 
MMHBOBJS=$(ZH2MOBJS) $(DH2MOBJS) $(DM2HOBJS) $(ZM2HOBJS)


all: dhb2mm zhb2mm dmm2hb zmm2hb

dhb2mm: $(DH2MOBJS)
	$(MPF90) -o dhb2mm   $(DH2MOBJS) $(PSBLAS_LIB) $(LDLIBS)
dmm2hb: $(DM2HOBJS)
	$(MPF90) -o dmm2hb   $(DM2HOBJS) $(PSBLAS_LIB) $(LDLIBS)
zhb2mm: $(ZH2MOBJS)
	$(MPF90) -o zhb2mm   $(ZH2MOBJS) $(PSBLAS_LIB) $(LDLIBS)
zmm2hb: $(ZM2HOBJS)
	$(MPF90) -o zmm2hb   $(ZM2HOBJS) $(PSBLAS_LIB) $(LDLIBS)

.f90.o:
	$(MPF90) $(F90COPT) $(FINCLUDES) $(FDEFINES) -c $<

clean: 
	/bin/rm -f $(MMHBOBJS) dhb2mm zhb2mm dmm2hb zmm2hb

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)

