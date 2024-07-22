include ../Make.inc


#CC=mpicc
#FC=mpif90
#FCOPT=-O0 -march=native
#OFFLOAD=-fopenacc -foffload=nvptx-none="-march=sm_70"

LIBDIR=../lib
INCDIR=../include
MODDIR=../modules
IMPLDIR=./impl 

INCLUDES=-I$(LIBDIR) -I$(INCDIR) -I$(MODDIR) 
FINCLUDES=$(FMFLAG). $(FMFLAG)$(INCDIR) $(FMFLAG)$(MODDIR) $(FIFLAG). 
CINCLUDES=
#LIBS=-L$(LIBDIR) -lpsb_util -lpsb_ext -lpsb_base -lopenblas -lmetis


FOBJS= psb_i_oacc_vect_mod.o \
       psb_s_oacc_vect_mod.o psb_s_oacc_csr_mat_mod.o \
       psb_d_oacc_vect_mod.o psb_d_oacc_csr_mat_mod.o \
       psb_c_oacc_vect_mod.o psb_c_oacc_csr_mat_mod.o \
       psb_z_oacc_vect_mod.o psb_z_oacc_csr_mat_mod.o \
	   psb_d_oacc_ell_mat_mod.o \
       psb_oacc_mod.o   psb_oacc_env_mod.o


LIBNAME=libpsb_openacc.a

OBJS=$(COBJS) $(FOBJS)


lib: objs ilib 
	ar cur $(LIBNAME) $(OBJS)
	/bin/cp -p $(LIBNAME) $(LIBDIR)

objs: $(OBJS) iobjs
	/bin/cp -p *$(.mod) $(MODDIR)

iobjs: $(OBJS)
	$(MAKE) -C impl objs 

ilib: $(OBJS)
	$(MAKE) -C impl lib 

psb_oacc_mod.o  : psb_i_oacc_vect_mod.o \
       psb_s_oacc_vect_mod.o psb_s_oacc_csr_mat_mod.o \
       psb_d_oacc_vect_mod.o psb_d_oacc_csr_mat_mod.o \
       psb_c_oacc_vect_mod.o psb_c_oacc_csr_mat_mod.o \
       psb_z_oacc_vect_mod.o psb_z_oacc_csr_mat_mod.o \
	   psb_d_oacc_ell_mat_mod.o \
       psb_oacc_env_mod.o

psb_s_oacc_vect_mod.o  psb_d_oacc_vect_mod.o \
       psb_c_oacc_vect_mod.o  psb_z_oacc_vect_mod.o : psb_i_oacc_vect_mod.o


psb_s_oacc_csr_mat_mod.o:        psb_s_oacc_vect_mod.o
psb_d_oacc_csr_mat_mod.o:        psb_d_oacc_vect_mod.o
psb_c_oacc_csr_mat_mod.o:        psb_c_oacc_vect_mod.o
psb_z_oacc_csr_mat_mod.o:        psb_z_oacc_vect_mod.o
psb_d_oacc_ell_mat_mod.o:        psb_d_oacc_vect_mod.o


clean: cclean iclean
	/bin/rm -f  $(FOBJS) *$(.mod) *.a
veryclean: clean
cclean: 
	/bin/rm -f  $(COBJS) 
iclean:
	$(MAKE) -C impl clean

.c.o:
	$(CC) $(CCOPT) $(CCOPENACC) $(CINCLUDES) $(CDEFINES) -c $< -o $@
.f90.o:
	$(FC) $(FCOPT) $(FCOPENACC) $(FINCLUDES) -c $< -o $@
.F90.o:
	$(FC) $(FCOPT) $(FCOPENACC) $(FINCLUDES) $(FDEFINES) -c $< -o $@
.cpp.o:
	$(CXX) $(CXXOPT) $(CXXOPENACC) $(CXXINCLUDES) $(CXXDEFINES) -c $< -o $@
