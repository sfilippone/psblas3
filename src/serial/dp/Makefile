include ../../../Make.inc
#
# The object files
#

FOBJS = dcrcr.o dcrdi.o dcrel.o dcrjd.o dgblock.o partition.o  \
	dgindex.o djadrp.o djadrp1.o dcsrrp.o dcsrp1.o check_dim.o \
	Max_nnzero.o dcoco.o dcocr.o dcrco.o djdcox.o djdco.o dvtfg.o dgind_tri.o \
	gen_block.o reordvn.o

# 
# dgind_tri.o
# 


OBJS=$(FOBJS)

#
# Where the library should go, and how it is called. 
# Note that we are regenerating most of libsparker.a on the fly. 
LIBDIR=../../../lib
#LIBNAME=libsparker.a	
LIBFILE=$(LIBDIR)/$(LIBNAME)
INCDIRS=-I. -I$(LIBDIR)

#
# No change should be needed below 
#


lib: $(FOBJS)

clean: cleanobjs

veryclean: cleanobjs

cleanobjs:
	/bin/rm -f $(OBJS) 


