.mod=.mod
.fh=.fh
.SUFFIXES: .f90 $(.mod)


#######################  Section 1 #######################
#  Define your compilers and compiler flags here         #
##########################################################
F90=/usr/local/gfortran/bin/gfortran
FC=/usr/local/gfortran/bin/gfortran
F77=$(FC)
CC=/usr/local/gfortran/bin/gcc
F90COPT= -O3  -ffast-math -march=pentium4 -msse2 -mfpmath=sse 
FCOPT=-O3  -ffast-math -march=pentium4 -msse2 -mfpmath=sse 
CCOPT=-O3 -ffast-math -march=pentium4 -msse2 -mfpmath=sse 

#######################  Section 2 #######################
#  Define your linker and linker flags here              #
##########################################################
F90LINK=/usr/local/mpich-gfortran/bin/mpif90 
FLINK=/usr/local/mpich-gfortran/mpif77 
MPF90=/usr/local/mpich-gfortran/bin/mpif90 
MPCC=/usr/local/mpich-gfortran/bin/mpicc 

#######################  Section 3 #######################
#  Specify paths to libraries                            #
##########################################################
BLAS=-lblas-gcc41 -L$(HOME)/LIB
BLACS=-lmpiblacs-gfortran -L$(HOME)/LIB


#######################  Section 4 #######################
#  Other useful tools&defines                            #
##########################################################
SLUDIR=/usr/local/SuperLU_3.0
SLU=-lslu_lx_gfort -L$(SLUDIR)
SLUDEF=-DHave_SLU_ -I$(SLUDIR)

UMFDIR=$(HOME)/LIB/Umfpack_gcc41
UMF=-lumfpack -lamd -L$(UMFDIR)
UMFDEF=-DHave_UMF_ -I$(UMFDIR)

#  Add -DLargeFptr for 64-bit addresses
CDEFINES=-DAdd_ $(SLUDEF) $(UMFDEF)

AR=ar -cur
RANLIB=ranlib



#######################  Section 5 #######################
#  Do not edit this                                      #
##########################################################
LIBDIR = lib
LIBNAME = libpsblas.a

TYPEMODS =  psb_spmat_type$(.mod) psb_descriptor_type$(.mod) psb_prec_type$(.mod) psb_realloc_mod$(.mod)
CONSTMODS = psb_const_mod$(.mod)
BLASMODS  = $(TYPEMODS) psb_psblas_mod$(.mod) psb_comm_mod$(.mod)
METHDMODS = psb_methd_mod$(.mod)
TOOLSMODS = $(CONSTMODS) psi_mod$(.mod) psb_tools_mod$(.mod) psb_serial_mod$(.mod)
PRECMODS  = psb_prec_mod$(.mod)
ERRORMODS = psb_error_mod$(.mod)
F90MODS= $(BLASMODS) $(PRECMODS) $(METHDMODS) $(TOOLSMODS) $(ERRORMODS) string$(.mod)

MODS=$(LIBDIR)/psb_const_mod$(.mod) $(LIBDIR)/psb_spmat_type$(.mod) $(LIBDIR)/psb_realloc_mod$(.mod) \
  $(LIBDIR)/psb_descriptor_type$(.mod)  $(LIBDIR)/psb_prec_type$(.mod) $(LIBDIR)/parts.fh \
  $(LIBDIR)/psb_serial_mod$(.mod) $(LIBDIR)/psb_comm_mod$(.mod) $(LIBDIR)/psb_error_mod$(.mod)

# Under Linux/gmake there is a rule interpreting .mod as Modula source! 
$(.mod).o:
	
.f.o:
	$(FC) $(FCOPT) $(INCDIRS)  -c $<
.c.o:
	$(CC) $(CCOPT) $(INCDIRS) $(CDEFINES) -c $<
.f$(.mod):
	$(F90) $(FCOPT) $(INCDIRS) -c $<
.f90$(.mod):
	$(F90) $(F90COPT) $(INCDIRS) -c $<
.f90.o:
	$(F90) $(F90COPT) $(INCDIRS) -c $<




