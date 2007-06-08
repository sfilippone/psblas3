.mod=.mod
.fh=.fh
.SUFFIXES: .f90 $(.mod) .F90 .F


#######################  Section 1 #######################
#  Define your compilers and compiler flags here         #
##########################################################
F90=/usr/local/gcc42/bin/gfortran
FC=/usr/local/gcc42/bin/gfortran
F77=$(FC)
CC=/usr/local/gcc42/bin/gcc
F90COPT=-O3 -pg
FCOPT=-O3   -pg
CCOPT=-O3   -pg
#
# Which flag does your compiler use  for module inclusion?
# Most compilers use -I  but Sun uses -M 
# FIFLAG is for INCLUDE 
FMFLAG=-I
FIFLAG=-I

#######################  Section 2 #######################
#  Define your linker and linker flags here              #
#  Does your MPI provide mpi.mod or mpif.h?              #
#   MPI_FINC=-DMPI_H  or MPI_FINC=-DMPI_MOD              #
#  If necessary add an include dir                       #
##########################################################
F90LINK=/usr/local/mpich-gcc42/bin/mpif90 -pg
FLINK=/usr/local/mpich-gcc42/bin/mpif77 
MPF90=/usr/local/mpich-gcc42/bin/mpif90 
MPF77=/usr/local/mpich-gcc42/bin/mpif77 	
MPCC=/usr/local/mpich-gcc42/bin/mpicc 
MPI_FINC=-DMPI_MOD
#$(FIFLAG)/opt/SUNWhpc/include

#######################  Section 3 #######################
#  Specify paths to libraries                            #
##########################################################
BLAS=-lblas-gcc42 -L$(HOME)/LIB
BLACS=-lmpiblacs-gcc42 -L$(HOME)/LIB


#######################  Section 4 #######################
#  Other useful tools&defines                            #
##########################################################
SLUDIR=/usr/local/SuperLU_3.0
SLU=-lslu_lx_gcc42 -L$(SLUDIR)
SLUDistDIR=/usr/local/SLUDist_2.0
SLUDist=-lslud-bck-gcc42 -L$(SLUDistDIR) -L $(HOME)/LIB
SLUDEF=-DHave_SLU_ -I$(SLUDIR) -DHave_SLUDist_ -I$(SLUDistDIR)

UMFDIR=$(HOME)/LIB/Umfpack_gcc41
UMF=-lumfpack -lamd -L$(UMFDIR)
UMFDEF=-DHave_UMF_ -I$(UMFDIR)
#
#  We are using the public domain tool METIS from U. Minnesota. To get it
#  check URL http://www.cs.umn.edu:~karypis
#  
METIS_LIB = -L$(HOME)/NUMERICAL/metis-4.0 -lmetis 
LDLIBS=$(BLACS) $(SLU) $(SLUDist) $(UMF) $(BLAS) $(METIS_LIB)

#  Add -DLargeFptr for 64-bit addresses
CDEFINES=-DAdd_ $(SLUDEF) $(UMFDEF)
FDEFINES=-DNETLIB_BLACS  $(MPI_FINC) -DHAVE_MOVE_ALLOC

AR=ar -cur
RANLIB=ranlib



#######################  Section 5 #######################
#  Do not edit this                                      #
##########################################################
LIBDIR=lib
BASELIBNAME=libpsb_base.a
PRECLIBNAME=libpsb_prec.a
METHDLIBNAME=libpsb_krylov.a
UTILLIBNAME=libpsb_util.a

# Under Linux/gmake there is a rule interpreting .mod as Modula source! 
$(.mod).o:
	
.c.o:
	$(CC) $(CCOPT) $(CINCLUDES) $(CDEFINES) -c $<
.f.o:
	$(FC) $(FCOPT) $(FINCLUDES)  -c $<
.f$(.mod):
	$(F90) $(FCOPT) $(FINCLUDES) -c $<
.f90$(.mod):
	$(F90) $(F90COPT) $(FINCLUDES) -c $<
.f90.o:
	$(F90) $(F90COPT) $(FINCLUDES) -c $<
.F.o:
	$(F90) $(F90COPT) $(FINCLUDES) $(FDEFINES) -c $<
.F90.o:
	$(F90) $(F90COPT) $(FINCLUDES) $(FDEFINES) -c $<
.F90$(.mod):
	$(F90) $(F90COPT) $(FINCLUDES) $(FDEFINES) -c $<




