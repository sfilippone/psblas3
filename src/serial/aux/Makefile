include ../../../Make.inc
#
# The object files
#

FOBJS  = isr.o isrx.o lsame.o \
	mrgsrt.o isaperm.o ibsrch.o


OBJS=$(FOBJS)

#
# Where the library should go, and how it is called. 
# Note that we are regenerating most of libsparker.a on the fly. 
LIBDIR=../../../lib
#LIBNAME=libsparker.a	
INCDIRS=-I. -I$(LIBDIR)

#
# No change should be needed below 
#

default: lib

lib: $(OBJS)
	$(AR) $(LIBDIR)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)

clean: 
	/bin/rm -f $(OBJS)

veryclean: clean

