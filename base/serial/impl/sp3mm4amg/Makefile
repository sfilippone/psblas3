include ../../../../Make.inc

CC=mpicc
CWALL =-Wall -Wextra
#extra to reduce useless warnings
CWALL+=-Wno-pointer-sign -Wno-unused-parameter -Wno-unused-but-set-variable -Wno-switch
CWALL+=-Wno-unused-function -Wno-unused-variable       #TODO CLEAN USELESS FUNCTIONS / VARIABLES
CWALL+=-Wno-unused-label -Wfatal-errors

CINCL   = -Iinclude/ 
CFLAGS  = -g -O3 $(CWALL) $(CINCL) -fopenmp $(RUNTIME)
MACROS = -DDEBUGPRINT="if(FALSE)" -DDEBUG="if(FALSE)" -DCONSISTENCY_CHECKS="if(FALSE)" -DVERBOSE="if(FALSE)"  -DDEBUGCHECKS="if(FALSE)" 
MACROSDBG = -DCONSISTENCY_CHECKS="if(TRUE)"  -DDEBUGCHECKS="if(TRUE)" -DVERBOSE="if(TRUE)" -DDEBUG="if(TRUE)" -DDEBUGPRINT="if(TRUE)"
LDFLAGS = -lm

LIBDIR=../../../
LIBFILE=$(LIBDIR)/$(LIBNAME)

SP3MM_CORE_OBJS=Sp3MM_CSR_OMP_Multi.o \
				Sp3MM_CSR_OMP_SymbStep_Multi.o \
				commons/ompGetICV.o \
				commons/sparseUtilsMulti.o \
				commons/utils.o \
				commons/ompi_dims_create/ompi_dims_create.o \
				lib/linuxK_rbtree_minimalized.o \
				lib/mmio.o \
				lib/parser.o

CBIND_OBJS=fbind/psb_f_spmm_ub.o

HEADERS = $(wildcard include/*.h)

OBJS= $(SP3MM_CORE_OBJS) $(CBIND_OBJS)

%.o : %.c $(HEADERS)
	$(CC) -c -o $@ $(CFLAGS) $< -Iinclude/ $(MACROSDBG) $(LDFLAGS)

objs: $(OBJS)

lib: objs
	$(AR) $(LIBDIR)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)

clean: cleanobjs

veryclean: cleanobjs

cleanobjs:
	/bin/rm -f $(OBJS) 
