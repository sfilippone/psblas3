INSTALLDIR=../../..
INCDIR=$(INSTALLDIR)/include/
MODDIR=$(INSTALLDIR)/modules/
include $(INCDIR)/Make.inc.psblas
#
# Libraries used
#
LIBDIR = $(INSTALLDIR)/lib/
PSBLAS_LIB = -L$(LIBDIR) -lpsb_util -lpsb_base
LDLIBS = $(PSBLDLIBS)

FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG).

EXEDIR=./runs


GREEN=\033[0;32m
RED=\033[0;31m
BLUE=\033[0;34m
YELLOW=\033[33m 
END_COLOUR=\033[0m


all: runsd psb_spmm_test
	@printf "$(GREEN)[INFO]\tCompilation success!$(END_COLOUR)\n" 

runsd: 	
	@(if test ! -d runs ; then mkdir runs; fi)
	@printf "$(BLUE)[INFO]\tBuild directory $(EXEDIR) correctly initialized$(END_COLOUR)\n"


psb_spmm_test:
	@$(FLINK) $(LOPT) psb_spmm_test.f90 spmm.f90 -o psb_spmm_test -I$(MODDIR) -I. $(PSBLAS_LIB) $(LDLIBS) 
	@mv psb_spmm_test $(EXEDIR)
	@printf "$(BLUE)[INFO]\tTesting files generated correctly$(END_COLOUR)\n" 

clean: 
	@rm -f $(OBJS)\
        *$(.mod) $(EXEDIR)/psb_spmm_test 

.PHONY: all runsd clean   

