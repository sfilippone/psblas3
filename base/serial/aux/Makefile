include ../../../Make.inc
#
# The object files
#

FOBJS  = isr.o isrx.o \
	mrgsrt.o isaperm.o ibsrch.o issrch.o  imsr.o imsrx.o


OBJS=$(FOBJS)

#
# Where the library should go, and how it is called. 
# Note that we are regenerating most of libsparker.a on the fly. 
SPARKERDIR=..
LIBDIR = ../..
MODDIR = ../../modules
INCDIRS = -I $(LIBDIR) -I $(MODDIR) -I . -I$(SPARKERDIR)

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

