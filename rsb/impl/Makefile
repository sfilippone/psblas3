include ../../Make.inc
PSBLIBDIR=$(PSBLASDIR)/lib/
PSBINCDIR=$(PSBLASDIR)/include
PSBMODDIR=$(PSBLASDIR)/modules
LIBDIR=../../lib
INCDIR=../../include
MODDIR=../../modules
PSBLAS_LIB= -L$(PSBLIBDIR)  -lpsb_util  -lpsb_base
#-lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base 
LDLIBS=$(PSBLDLIBS)
#
# Compilers and such
#
#CCOPT= -g
FINCLUDES=$(FMFLAG).. $(FMFLAG)$(INCDIR) $(FMFLAG)$(MODDIR) $(FMFLAG)$(PSBMODDIR) $(FMFLAG)$(PSBINCDIR) $(LIBRSB_INCLUDES) $(FIFLAG).. $(LIBRSB_DEFINES)
CINCLUDES=
LIBNAME=libpsb_rsb.a

OBJS= \
psb_d_cp_rsb_from_coo.o \
psb_d_mv_rsb_from_coo.o \
psb_d_cp_rsb_to_coo.o psb_d_rsb_csmv.o

objs: $(OBJS)

lib: objs
	ar cur ../$(LIBNAME) $(OBJS)

clean:
	/bin/rm -f $(OBJS)
