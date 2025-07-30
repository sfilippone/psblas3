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
EXEDIR = ./runs

FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG).

GREEN=\033[0;32m
RED=\033[0;31m
BLUE=\033[0;34m
YELLOW=\033[33m 
END_COLOUR=\033[0m


all: runsd psb_geaxpby_test 
	@printf "$(GREEN)[INFO]\t  Compilation success!$(END_COLOUR)\n" 

runsd: 	
	@(if test ! -d $(EXEDIR) ; then mkdir $(EXEDIR); fi)
	@printf "$(BLUE)[INFO]\t  Build directory $(EXEDIR) correctly initialized$(END_COLOUR)\n"

psb_geaxpby_test:
	@$(FLINK) $(LOPT) psb_geaxpby_test.f90 ../utils/psb_test_utils.o -o $(EXEDIR)/psb_geaxpby_test -I../utils/ -I$(MODDIR) -I. $(PSBLAS_LIB) $(LDLIBS) 
	@printf "$(BLUE)[INFO]\t  Testing files for psb_geaxpby linked and compiled correctly$(END_COLOUR)\n" 

clean: 
	@rm -f $(OBJS)\
        *$(.mod) $(EXEDIR)/psb_geaxpby_test 

.PHONY: all runsd clean   
