## $Id$
##---------------------------------------------------------------------------
## LaTeX Makefile
## Copyright (C) 1996-2001  Michael Forman	Michael.Forman@Colorado.EDU
## 
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
##
## This copyright applies to this Makefile, and all perl scripts.  
## The GPL does NOT apply to the actual content of the paper or thesis.  
##---------------------------------------------------------------------------
## 01.Dec,1996	forman	 Initial Makefile 
## 01.Jun,1997	forman	 Added support for print, preview, and bibtex
## 01.Jan,1998	stones	 tmp and lib directories to reduce clutter
## 05.Feb,1998	forman	 Added the following functions:
##			 .PHONY to prevent errors
##			 generic TEXFILE definition with following patsubst's
##			 vpath definitions
##			 gzip, gunzip, tar, ci
##			 search for \bibliography in tex file
##			 conditional ifeq omits bibtex if unnecessary
##			 documentation added
## 30.Nov,1998	forman	 Added TOPFILE and SECFILE definitions to support
##        		 texfiles with "input" commands.
## 23.Dec,1998	marshats Added RCSFILES and ci/co capabilites for multiple files
##                       Now only runs makeindex if $(IDX) file exists.
## 15.Mar,1999	forman	 Added 'make wc' for papers with word quotas.
## 24.Aug,1999	forman	 Converted Makefile to use pdftex as the primary
##			 complier.  Now generates true pdf and ps files.
## 09.Feb,2000	forman	 Simplified the documentation.  Removed direct
##			 compilation of tex into DVI and PS.  Conversion
##			 is done with PDF2PS now.
## 08.Aug,2000	forman	 Added define for figures directories.
## 18.Aug,2000	forman	 Removed the redundant code in the $(pdflatex) and
##			 $(pdflatex-bibtex) definitions by splitting them up
##			 into several smaller definitions which are called
##			 from a single set of "if-then" statements.
## 18.Aug,2000	forman	 Added support for glosstex and makeindex.
## 19.Aug,2000	forman	 To save space, all files in tmp are now links.
## 19.Aug,2000	forman	 Removed all postscript commands in favor of pdf.
##---------------------------------------------------------------------------
##
## This Makefile expects the following directory structure:
##
##     	Makefile	This file
##     	*.tex		Put latex files in root directory.
##     	RCS/		Create an RCS directory for "ci" and "co".
##	doc/		Document directory.  Includes gpl.txt.
##     	lib/		Put all cls, sty, idx, gdf, and bib files in lib.
##     	figures/	Put all figures in the figures directory.
##     	tmp/		Never put anything in tmp -- it gets cleaned out.
##
##---------------------------------------------------------------------------
##
## Normal Usage:
##	make 		Run pdflatex
##	make clean 	Remove all files in tmp and the pdf link in root.
##	make preview 	Preview the compiled file
##	make ci 	Check in the RCSFILES
##	make co 	Check out the RCSFILES
##
## Advanced Usage:
##	make gzip	Recursively gzip all the files in the root directory
##	make gunzip	Reverse the above process
##	make tar 	Tar and gzip the working directory
##	make wc		Count the number of words in your report
##
##---------------------------------------------------------------------------
## WARNING:
## 	If "make ci" fails, "make co" will nuke your files!
## 	Don't put anything in tmp, "make clean" will delete everything!
##---------------------------------------------------------------------------
#

TOPFILE   = userguide.tex
SECFILE   = title.tex intro.tex commrout.tex datastruct.tex psbrout.tex toolsrout.tex\
	methods.tex precs.tex penv.tex error.tex util.tex
FIGDIR    = figures

XPDFFLAGS = 
ACROFLAGS = 				#-- macos
#ACROFLAGS = -geometry 1234x1168+0+0	#-- new-sydney-wide
#ACROFLAGS = -geometry 1000x1000+0+0	#-- sydney-wide
#ACROFLAGS = -geometry 750x1000+0+0	#-- sydney
#ACROFLAGS = -geometry 1200x1200+0+0	#-- home-wide
#ACROFLAGS = -geometry 900x1200+0+0	#-- home

##---------------------------------------------------------------------------
## Change nothing below here (unless you're really really good).
#----------------------------------------------------------------------------

##
## Accounting
#
START    = $(shell date)
WDIR     = $(notdir $(shell pwd))
DATE     = $(shell date +%Y-%m-%d)

##
## Programs
#
PDFLATEX = pdflatex
ACRO     = evince
XPDF     = xpdf
WC       = wc
PDF2PS   = pdf2ps
PDF2TEXT = pdftotext
MAKEIDX  = makeindex
GLOSSTEX = glosstex
BIBTEX   = bibtex
FILTER   = ../bin/texfilter 
FILTER   = cat
CLEANIDX = ../bin/cleanidx

##
## Files
#
TEXFILES = $(TOPFILE) $(SECFILE)
RCSFILES = $(TEXFILES) Makefile

BASEFILE = $(patsubst %.tex,%,$(TOPFILE))
IDX      = $(join $(BASEFILE),.idx)
PDF      = $(join $(BASEFILE),.pdf)
PS       = $(join $(BASEFILE),.ps)
GXS      = $(join $(BASEFILE),.gxs)
GLX      = $(join $(BASEFILE),.glx)

FIGURES  = $(sort $(wildcard $(FIGDIR)/*))
GLOFILES:= $(sort $(wildcard lib/*.gdf))
GLOFILES:= $(patsubst lib/%,%,$(GLOFILES))
BIBFILES:= $(shell grep ^.bibliography{ $(TOPFILE)|sed "s/,/ /g"|sed "s/%.*//g")
BIBFILES:= $(patsubst \bibliography{%,%,$(BIBFILES))
BIBFILES:= $(patsubst %},%,$(BIBFILES))
BIBFILES:= $(sort $(addsuffix .bib,$(BIBFILES)))

LIBFILES = $(sort $(wildcard lib/*))

TEXLNFIL = $(addprefix $(PWD)/,$(TEXFILES))

#============================================================================

all: $(PDF)
	@echo $(PDF)

$(PDF): $(TEXFILES) $(LIBFILES) $(FIGURES) Makefile
	$(header)
	$(initialize)
	$(pdflatex-filter)
        ifneq ($(BIBFILES),)
	  $(bibtex)
	  $(pdflatex-filter)
        endif
        ifneq ($(GLOFILES),)
	  $(glosstex)
	  $(pdflatex-filter)
        endif
#	$(makeindex)
	$(pdflatex-filter)
	$(finish)

$(PS): $(PDF)
	$(PDF2PS) $< $(PS)

ps: $(PS) 

#============================================================================

.PHONY: clean clean-tmp preview print gzip gunzip tar ci

clean:
	rm -f $(PDF) $(PS)
	cd tmp ; rm -f *

clean-tmp:
	cd tmp ; rm -f *

#----------------------------------------------------------------------------

preview: $(PDF)
	$(ACRO) $(ACROFLAGS) $< &

xpreview: $(PDF)
	$(XPDF) $(XPDFFLAGS) $< &

wc: $(PDF)
	$(PDF2TEXT) $< | $(WC) 

#----------------------------------------------------------------------------

ci: $(RCSFILES)
	$(shell for i in $(RCSFILES) ; do ci -u $$i ; done) 

co: $(RCSFILES)
	$(shell for i in $(RCSFILES) ; do co -l $$i ; done) 

gzip: clean-tmp
	gzip -r ./
	gunzip Makefile.gz

gunzip: 
	gunzip -r ./

tar: clean-tmp
	cd ..                               ; \
	tar cvf $(WDIR)-$(DATE).tar $(WDIR) ; \
	gzip $(WDIR)-$(DATE).tar

info: 
	$(header)

#============================================================================

define header
  @echo
  @echo "#---------------------------------------------------------------------"
  @echo "MAKEFILE = LaTeX PDF Makefile"
  @echo "AUTHOR   = Alfredo Buttari"
  @echo 'ID       = $$Id$ ' 
  @echo "#---------------------------------------------------------------------"
  @echo
  @echo "ACRO     = $(ACRO) $(ACROFLAGS) $(PDF)"
  @echo "XPDF     = $(XPDF) $(XPDFFLAGS) $(PDF)"
  @echo "GV       = $(GV) $(GVFLAGS) $(PS)"
  @echo "LPR      = $(LPR) $(LPRFLAGS) $(PS)"
  @echo
  @echo "WDIR     = $(WDIR)"
  @echo "DATE     = $(DATE)"
  @echo
  @echo "TOPFILE  = $(TOPFILE)"
  @echo "SECFILE  = $(SECFILE)"
  @echo "TEXFILES = $(TEXFILES)"
  @echo "PDF      = $(PDF)"
  @echo "PS       = $(PS)"
  @echo "BIBFILES = $(BIBFILES)"
  @echo "GLOFILES = $(GLOFILES)"
  @echo "IDXFILES = $(IDXFILES)"
  @echo
endef
# @echo "FIGURES  = $(FIGURES)"

#----------------------------------------------------------------------------


define initialize
  @if test ! -d tmp; then mkdir tmp; fi
  @ln -sf $(TEXLNFIL) tmp
  @ln -sf $(PWD)/lib/* tmp
  @ln -sf $(PWD)/$(FIGDIR) tmp
endef

#----------------------------------------------------------------------------

define pdflatex
  @echo
  @echo "----- pdflatex -------------------------------------------------------"
  @echo -n "Starting: "; date
  @echo
  cd tmp; $(PDFLATEX) $(TOPFILE)
endef

#----------------------------------------------------------------------------

define pdflatex-filter
  @echo
  @echo "----- pdflatex -------------------------------------------------------"
  @echo -n "Starting: "; date
  @echo
  cd tmp; ($(PDFLATEX) $(TOPFILE) 2>&1) | $(FILTER)
endef

#----------------------------------------------------------------------------

define bibtex
  @echo
  @echo "----- bibtex ---------------------------------------------------------"
  @echo -n "Starting: "; date
  @echo
  cd tmp; $(BIBTEX) $(BASEFILE)
endef

#----------------------------------------------------------------------------

define glosstex
  @echo
  @echo "----- glosstex -------------------------------------------------------"
  @echo -n "Starting: "; date
  @echo
  cd tmp; $(GLOSSTEX) $(BASEFILE) $(GLOFILES)
  cd tmp; $(MAKEIDX) $(GXS) -o $(GLX) -s glosstex.ist
endef

#----------------------------------------------------------------------------

define makeindex
  @echo
  @echo "----- makeindex ------------------------------------------------------"
  @echo -n "Starting: "; date
  @echo
  cd tmp; mv $(IDX) $(IDX)-; $(CLEANIDX) < $(IDX)- > $(IDX)
  cd tmp; $(MAKEIDX) $(IDX)
endef

#----------------------------------------------------------------------------

define finish
  @ln -sf tmp/$@ .
  @echo
  @echo "----- finish ---------------------------------------------------------"
  @echo -n "Start:  "$(START); echo
  @echo -n "Finish: "; date
  @echo -n "Output: "; ls -l -o tmp/$@
  @echo -n "Target: "
endef


