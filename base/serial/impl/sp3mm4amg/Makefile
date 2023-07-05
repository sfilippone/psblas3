include ../../../../Make.inc

CC=gcc
CWALL =-Wall -Wextra
#extra to reduce useless warnings
CWALL+=-Wno-pointer-sign -Wno-unused-parameter -Wno-unused-but-set-variable -Wno-switch
CWALL+=-Wno-unused-function -Wno-unused-variable       #TODO CLEAN USELESS FUNCTIONS / VARIABLES
CWALL+=-Wno-unused-label -Wfatal-errors

CINCL   = -Iinclude/ 
CFLAGS  = -O2 $(CWALL) $(CINCL) -fopenmp $(RUNTIME)
LDFLAGS = -lm

LIBDIR=../../../
LIBFILE=$(LIBDIR)/$(LIBNAME)

SP3MM_CORE_OBJS=commons/ompGetICV.o \
				commons/sparseUtilsMulti.o \
				commons/utils.o \
				lib/linuxK_rbtree_minimalized.o \
				lib/mmio.o \
				lib/parser.o

CBIND_OBJS=fbind/psb_f_spmm_ub.o

OBJS= $(SP3MM_CORE_OBJS) $(CBIND_OBJS)

%.o : %.c
	$(CC) -c -o $@ $(CFLAGS) $^ $(LDFLAGS)

objs: $(OBJS)

lib: objs
	$(AR) $(LIBDIR)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)

clean: cleanobjs

veryclean: cleanobjs

cleanobjs:
	/bin/rm -f $(OBJS) 
