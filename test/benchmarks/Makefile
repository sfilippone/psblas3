INSTALLDIR=../..
INCDIR=$(INSTALLDIR)/include
MODDIR=$(INSTALLDIR)/modules
include $(INCDIR)/Make.inc.psblas
#
# Libraries used
#
LIBDIR=$(INSTALLDIR)/lib
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base
LDLIBS=$(PSBLDLIBS)

CCOPT= -g
FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG).

SFOBJS=getp.o  psb_sf_sample.o
DFOBJS=getp.o  psb_df_sample.o
CFOBJS=getp.o  psb_cf_sample.o
ZFOBJS=getp.o  psb_zf_sample.o
DTOBJS=d_file_spmv.o
STOBJS=s_file_spmv.o
DPGOBJS=pdgenspmv.o

EXEDIRFILEREAD=./runs-fileread
EXEDIRKERNEL=./runs-kernel
EXEDIRPARGEN=./runs-pargen

all:    psb_sf_sample psb_df_sample psb_cf_sample psb_zf_sample \
	d_file_spmv s_file_spmv pdgenspmv \
	psb_d_pde3d psb_s_pde3d psb_d_pde2d psb_s_pde2d

psb_sf_sample.o psb_df_sample.o psb_cf_sample.o psb_zf_sample.o: getp.o

psb_sf_sample: $(SFOBJS)
	$(FLINK) $(LOPT) $(SFOBJS) -o psb_sf_sample $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_sf_sample $(EXEDIR-FILEREAD)

psb_df_sample: $(DFOBJS)
	$(FLINK) $(LOPT) $(DFOBJS) -o psb_df_sample $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_df_sample $(EXEDIR-FILEREAD)

psb_cf_sample: $(CFOBJS)
	$(FLINK) $(LOPT) $(CFOBJS) -o psb_cf_sample $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_cf_sample $(EXEDIR-FILEREAD)

psb_zf_sample: $(ZFOBJS)
	$(FLINK) $(LOPT) $(ZFOBJS) -o psb_zf_sample $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_zf_sample $(EXEDIR-FILEREAD)

d_file_spmv: $(DTOBJS)
	$(FLINK) $(LOPT) $(DTOBJS) -o d_file_spmv $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv d_file_spmv $(EXEDIR-KERNEL)

pdgenspmv: $(DPGOBJS)
	$(FLINK) $(LOPT) $(DPGOBJS) -o pdgenspmv $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv pdgenspmv $(EXEDIR-KERNEL)

s_file_spmv: $(STOBJS)
	$(FLINK) $(LOPT) $(STOBJS) -o s_file_spmv $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv s_file_spmv $(EXEDIR-KERNEL)

psb_d_pde3d: psb_d_pde3d.o
	$(FLINK)  psb_d_pde3d.o -o psb_d_pde3d $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_d_pde3d $(EXEDIR-PARGEN)

psb_s_pde3d: psb_s_pde3d.o
	$(FLINK) psb_s_pde3d.o  -o psb_s_pde3d $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_s_pde3d $(EXEDIR-PARGEN)

psb_d_pde2d: psb_d_pde2d.o
	$(FLINK)  psb_d_pde2d.o -o psb_d_pde2d $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_d_pde2d $(EXEDIR-PARGEN)

psb_s_pde2d: psb_s_pde2d.o
	$(FLINK) psb_s_pde2d.o  -o psb_s_pde2d $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_s_pde2d $(EXEDIR-PARGEN)


.f90.o:
	$(MPFC) $(FCOPT) $(FINCLUDES) $(FDEFINES) -c $<

clean:
	/bin/rm -f $(DFOBJS) $(ZFOBJS) $(SFOBJS) $(CFOBJS) \
	$(DBOBJSS) $(DBOBJS) $(DTOBJS) $(STOBJS) \
	psb_d_pde3d.o psb_s_pde3d.o  psb_d_pde2d.o psb_s_pde2d.o  \
	$(EXEDIR-PARGEN)/psb_d_pde3d $(EXEDIR-PARGEN)/psb_s_pde3d \
	$(EXEDIR-PARGEN)/psb_d_pde2d $(EXEDIR-PARGEN)/psb_s_pde2d \
        *$(.mod) $(EXEDIR-FILEREAD)/psb_*f_sample

lib:
	(cd ../../; make library)

verycleanlib:
	(cd ../../; make veryclean)
