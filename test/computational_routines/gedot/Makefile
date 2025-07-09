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

GREEN=\033[0;32m
RED=\033[0;31m
BLUE=\033[0;34m
YELLOW=\033[33m 
END_COLOUR=\033[0m


all: psb_gedot_test 
	@printf "$(GREEN)[INFO]\t  Compilation success!$(END_COLOUR)\n" 

psb_gedot_test:
	@$(FLINK) $(LOPT) psb_gedot_test.f90 ../utils/psb_test_utils.o -o psb_gedot_test -I../utils/ -I$(MODDIR) -I. $(PSBLAS_LIB) $(LDLIBS) 
	@printf "$(BLUE)[INFO]\t  Testing files for psb_gedot linked and compiled correctly$(END_COLOUR)\n" 


clean: 
	@rm -f $(OBJS)\
        *$(.mod) psb_gedot_test 

.PHONY: all psb_gedot_test clean   
