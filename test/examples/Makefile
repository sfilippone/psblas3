BASEDIR=../..
INCDIR=$(BASEDIR)/include
MODDIR=$(INSTALLDIR)/modules/
include $(INCDIR)/Make.inc.psblas
#
# Libraries used
#
LIBDIR=$(BASEDIR)/lib
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base
LDLIBS=$(PSBLDLIBS)
#
# Compilers and such
#
CCOPT= -g
FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG).
EXEDIR=./runs


all: hello pingpong d_matgen

hello: hello.o
	$(FLINK)  hello.o -o hello $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv hello $(EXEDIR)

pingpong: pingpong.o
	$(FLINK)  pingpong.o -o pingpong $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv pingpong $(EXEDIR)

psb_d_xyz_impl.o d_matgen.o: psb_d_xyz_mat_mod.o

d_matgen: d_matgen.o psb_d_xyz_mat_mod.o psb_d_xyz_impl.o
	$(FLINK) $(LINKOPT) d_matgen.o psb_d_xyz_mat_mod.o psb_d_xyz_impl.o \
	 -o d_matgen  $(PSBLAS_LIB) $(LDLIBS)
	/bin/cp -p $(CPUPDFLAG) d_matgen $(EXEDIR)
#	/bin/mv d_matgen $(EXEDIR)

check: all
	cd runs && echo 5 | ./d_matgen

clean:
	/bin/rm -f hello.o pingpong.o \
	d_matgen.o psb_d_xyz_mat_mod.o psb_d_xyz_impl.o \
	$(EXEDIR)/hello *$(.mod)

verycleanlib:
	(cd ../..; make veryclean)

lib:
	(cd ../../; make library)
