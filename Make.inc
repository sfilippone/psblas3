.mod=.mod
.fh=.fh
.SUFFIXES: .f90 $(.mod)


#######################  Section 1 #######################
#  Define your compilers and compiler flags here         #
##########################################################
F90=ifort
FC=ifort
CC=icc
F77=$(FC)
F90COPT= -check all -g -CB -no_cpprt -warn interfaces
FCOPT=-check all -g -CB -no_cpprt -warn interfaces
CCOPT=

#######################  Section 2 #######################
#  Define your linker and linker flags here              #
##########################################################
F90LINK=/usr/local/mpich-intel/bin/mpif90 -g -CB -no_cpprt
FLINK=mpif77 -g -CB -no_cpprt
MPF90=/usr/local/mpich-intel/bin/mpif90 -g -CB -no_cpprt
MPCC=/usr/local/mpich-intel/bin/mpicc -g -CB -no_cpprt

#######################  Section 3 #######################
#  Specify paths to libraries                            #
##########################################################
BLAS=-lblas-intel -L$(HOME)/NUMERICAL/LIB
BLACS=-lmpiblacs-intel -L$(HOME)/NUMERICAL/LIB


#######################  Section 4 #######################
#  Other useful tools&defines                            #
##########################################################
CDEFINES=-DAdd_
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




