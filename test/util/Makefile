BASEDIR=../..
INCDIR=$(BASEDIR)/include/
MODDIR=$(INSTALLDIR)/modules/
include $(INCDIR)/Make.inc.psblas
#
# Libraries used
#
LIBDIR=$(BASEDIR)/lib/
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base
LDLIBS=$(PSBLDLIBS)

FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG).

ZH2MOBJS=zhb2mm.o 
DH2MOBJS=dhb2mm.o 
DM2HOBJS=dmm2hb.o 
ZM2HOBJS=zmm2hb.o 
MMHBOBJS=$(ZH2MOBJS) $(DH2MOBJS) $(DM2HOBJS) $(ZM2HOBJS)


all: dhb2mm zhb2mm dmm2hb zmm2hb

dhb2mm: $(DH2MOBJS)
	$(MPFC) -o dhb2mm   $(DH2MOBJS) $(PSBLAS_LIB) $(LDLIBS)
dmm2hb: $(DM2HOBJS)
	$(MPFC) -o dmm2hb   $(DM2HOBJS) $(PSBLAS_LIB) $(LDLIBS)
zhb2mm: $(ZH2MOBJS)
	$(MPFC) -o zhb2mm   $(ZH2MOBJS) $(PSBLAS_LIB) $(LDLIBS)
zmm2hb: $(ZM2HOBJS)
	$(MPFC) -o zmm2hb   $(ZM2HOBJS) $(PSBLAS_LIB) $(LDLIBS)

.f90.o:
	$(MPFC) $(FCOPT) $(FINCLUDES) $(FDEFINES) -c $<

clean: 
	/bin/rm -f $(MMHBOBJS) dhb2mm zhb2mm dmm2hb zmm2hb

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)

