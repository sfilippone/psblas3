include ../../../../Make.inc
#
# The object files
#

FOBJS = dcrcr.o dcrdi.o dcrel.o dcrjd.o dgblock.o partition.o  \
	dgindex.o djadrp.o djadrp1.o dcsrrp.o dcsrp1.o check_dim.o \
	Max_nnzero.o dcoco.o dcocr.o dcrco.o dcrinco.o djdcox.o djdco.o dvtfg.o dgind_tri.o \
	gen_block.o dcoinco.o reordvn.o zreordvn.o\
	zcrcr.o zcsrrp.o zcsrp1.o zgindex.o zgind_tri.o zcocr.o zcrinco.o zcoco.o

# 
# dgind_tri.o
# 


OBJS=$(FOBJS)

#
# Where the library should go, and how it is called. 
# Note that we are regenerating most of libsparker.a on the fly. 
#LIBDIR=../../LIB
#LIBNAME=libsparker.a	
LIBFILE=$(LIBDIR)/$(LIBNAME)
SPARKERDIR=..
INCDIRS=-I. -I$(SPARKERDIR) -I$(LIBDIR)

#
# No change should be needed below 
#


lib: $(FOBJS)
	$(AR) $(LIBFILE) $(OBJS) 
	$(RANLIB) $(LIBFILE)

$(FOBJS): $(SPARKERDIR)/sparker.fh

clean: cleanobjs

veryclean: cleanobjs

cleanobjs:
	/bin/rm -f $(OBJS) 


