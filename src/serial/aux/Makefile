include ../../../../Make.inc
#
# The object files
#

FOBJS  = daxpby.o    	getrepflag.o 	geterr.o   \
	isr.o		isrx.o		lsame.o    \
	setrepflag.o	seterr.o	sperror.o  \
	write_message.o	mrgsrt.o	xerbla.o   \
	xsperr.o	zaxpby.o	zseterr.o  \
	zsperror.o	zwrite_message.o                   \
	zxsperr.o    zsetrepflag.o  isaperm.o ibsrch.o


OBJS=$(FOBJS)

#
# Where the library should go, and how it is called. 
# Note that we are regenerating most of libsparker.a on the fly. 
#LIBDIR=../../../LIB
#LIBNAME=libsparker.a	
LIBFILE=$(LIBDIR)/$(LIBNAME)
SPARKERDIR=..
INCDIRS=-I. -I$(SPARKERDIR) -I$(LIBDIR)

#
# No change should be needed below 
#

default: lib

lib: $(OBJS)
	$(AR) $(LIBFILE) $(OBJS) 
	$(RANLIB) $(LIBFILE)

$(FOBJS): $(SPARKERDIR)/sparker.fh

clean: cleanobjs

veryclean: cleanobjs

cleanobjs:
	/bin/rm -f $(OBJS) 

