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


all: runsd psb_geaxpby_test 
	@printf "$(GREEN)[INFO]\t  Compilation success!$(END_COLOUR)\n" 

check_changes:
	@echo "Checking for changes..."
	@if [ "$(find . -name '*.f90' -o -name '*.h' -o -name '*.c' -newer build)" ]; then \
		echo "Changes detected. Running make..."; \
		$(MAKE) build; \
	else \
		echo "No changes detected."; \
	fi


runsd: 	
	@(if test ! -d runs ; then mkdir runs; fi)
	@printf "$(BLUE)[INFO]\t  Build directory $(EXEDIR) correctly initialized$(END_COLOUR)\n"


psb_geaxpby_test:
	@$(FLINK) $(LOPT) psb_geaxpby_test.f90 geaxpby.f90 -o psb_geaxpby_test -I$(MODDIR) -I. $(PSBLAS_LIB) $(LDLIBS) 
	@mv psb_geaxpby_test $(EXEDIR)
	@printf "$(BLUE)[INFO]\t  Testing files generated correctly$(END_COLOUR)\n" 


clean: 
	@rm -f $(OBJS)\
        *$(.mod) $(EXEDIR)/psb_geaxpby_test 

.PHONY: all runsd clean   
