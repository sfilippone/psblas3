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


all: 
	@$(FLINK) $(LOPT) -c psb_test_utils.f90 -o psb_test_utils.o -I$(MODDIR) -I. $(PSBLAS_LIB) $(LDLIBS) 
	@printf "$(BLUE)[INFO]\t  Utilities for testing files generated correctly$(END_COLOUR)\n" 
	@printf "$(GREEN)[INFO]\t  Compilation success!$(END_COLOUR)\n" 


clean: 
	@rm -f *.o *$(.mod)  

.PHONY: all clean   