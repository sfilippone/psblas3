#-----------------------------------
# Set oldest allowable CMake version
#-----------------------------------
cmake_minimum_required(VERSION 3.4)

#----------------------------------------------
# Define canonical CMake build types and extras
#----------------------------------------------
set ( CMAKE_CONFIGURATION_TYPES "Debug" "Release" "MinSizeRel" "RelWithDebInfo" "CodeCoverage" )
set ( CMAKE_BUILD_TYPE "Release"
  CACHE STRING "Select which configuration to build." )
set_property ( CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${CMAKE_CONFIGURATION_TYPES} )

#-----------------------------------------------------
# Determine version from .VERSION file or git describe
#-----------------------------------------------------
include("${CMAKE_SOURCE_DIR}/cmake/setVersion.cmake")
set_version(
  VERSION_VARIABLE PSBLAS_Version
  GIT_DESCRIBE_VAR full_git_describe
  CUSTOM_VERSION_FILE "${CMAKE_SOURCE_DIR}/.VERSION")
message( STATUS "Building PSBLAS version: ${full_git_describe}" )

#------------------------------------------
# Name project and specify source languages
#------------------------------------------
project(psblas VERSION "${PSBLAS_Version}" LANGUAGES C Fortran)

#-----------------------------------------------------------------
# Define a target to create a checksummed & signed release archive
#-----------------------------------------------------------------
set(${CMAKE_PROJECT_NAME}_dist_string "${CMAKE_PROJECT_NAME}-${full_git_describe}")
if(GIT_FOUND)
  add_custom_target(dist # OUTPUT "${CMAKE_BINARY_DIR}/${_package_stem_name}.tar.gz"
    COMMAND "${CMAKE_COMMAND}" -P "${CMAKE_SOURCE_DIR}/cmake/makeDist.cmake" "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}"
    COMMENT "Creating source release asset, ${_package_stem_name}.tar.gz, from ${_full_git_describe} using the `git archive` command."
    VERBATIM)
endif()

#--------------------------
# Prohibit in-source builds
#--------------------------
if ("${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_BINARY_DIR}")
  message(FATAL_ERROR "ERROR! "
    "CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}"
    " == CMAKE_CURRENT_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}"
    "\nThis archive does not support in-source builds:\n"
    "You must now delete the CMakeCache.txt file and the CMakeFiles/ directory under "
    "the 'src' source directory or you will not be able to configure correctly!"
    "\nYou must now run something like:\n"
    "  $ rm -r CMakeCache.txt CMakeFiles/"
    "\n"
    "Please create a directory outside the ${CMAKE_PROJECT_NAME} source tree and build under that outside directory "
    "in a manner such as\n"
    "  $ mkdir build\n"
    "  $ cd build\n"
    "  $ CC=gcc FC=gfortran cmake -DBUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/path/to/install/dir /path/to/psblas3/src/dir \n"
    "\nsubstituting the appropriate syntax for your shell (the above line assumes the bash shell)."
    )
endif()

#----------------------------------------------------
# Define coverage flags and report untested compilers
#----------------------------------------------------
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU" )
  set(gfortran_compiler true)
  set ( CMAKE_C_FLAGS_CODECOVERAGE "-fprofile-arcs -ftest-coverage -O0"
    CACHE STRING "Code coverage C compiler flags")
  set ( CMAKE_Fortran_FLAGS_CODECOVERAGE "-fprofile-arcs -ftest-coverage -O0"
    CACHE STRING "Code coverage Fortran compiler flags")
else()
  message(WARNING
    "\n"
    "Attempting untested CMake build with Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}. "
    "Please report any failures at https://github.com/sfilippone/psblas3\n\n"
  )
endif()

#-------------------------------------
# Fortran name mangling intro-spection
#-------------------------------------
include("${CMAKE_SOURCE_DIR}/cmake/CapitalizeString.cmake")
include(FortranCInterface)
CapitalizeString(${FortranCInterface_GLOBAL__CASE} fc_case)
message(STATUS "Name mangling capitalization: ${fc_case}")
message(STATUS "Name mangling fortran global suffix underscore: ${FortranCInterface_GLOBAL__SUFFIX}")
if(FortranCInterface_GLOBAL__SUFFIX STREQUAL "")
  add_definitions("-D${fc_case}Case")
elseif(FortranCInterface_GLOBAL__SUFFIX STREQUAL "_")
  add_definitions("-D${fc_case}Underscore")
elseif(FortranCInterface_GLOBAL__SUFFIX STREQUAL "__")
  add_definitions("-D${fc_case}DoubleUnderscore")
else()
  message( FATAL_ERROR "Fortran name mangling suffix, \'${FortranCInterface_GLOBAL__SUFFIX}\', unknown to morfeus/NEMO")
endif()

#----------------------------------------------
# Determine system endian-ness and pointer size
#----------------------------------------------
include(TestBigEndian)
TEST_BIG_ENDIAN(IS_BIG_ENDIAN)
if(IS_BIG_ENDIAN)
  message( STATUS "System appears to be big endian.")
else()
  message( STATUS "System appears to be little endian.")
  add_definitions(-DLittleEndian)
endif()
include(CheckTypeSize)
CHECK_TYPE_SIZE("void *" VOID_P_SIZE LANGUAGE C)
if(${VOID_P_SIZE} EQUAL 8)
  add_definitions(-DPtr64Bits)
endif()
message(STATUS "Have 64bit pointers")

#-----------------------------------------
# Check for some Fortran compiler features
#-----------------------------------------
include(CheckFortranSourceCompiles)
CHECK_Fortran_SOURCE_COMPILES(
  "
integer, allocatable :: a(:), b(:)
allocate(a(5))
a = [1,2,3,4,5]
call move_alloc(from=a, to=b)
end
"
  HAVE_MOVE_ALLOC
  SRC_EXT f90
  )
if(HAVE_MOVE_ALLOC)
  add_definitions(-DHAVE_MOVE_ALLOC)
endif()
CHECK_Fortran_SOURCE_COMPILES(
  "integer, volatile :: i ; end"
  HAVE_VOLATILE
  SRC_EXT f90
  )
if(HAVE_VOLATILE)
  add_definitions(-DHAVE_VOLATILE)
endif()
CHECK_Fortran_SOURCE_COMPILES(
  "use ISO_FORTRAN_ENV ; end"
  HAVE_ISO_FORTRAN_ENV
  SRC_EXT f90
  )
if(HAVE_ISO_FORTRAN_ENV)
  add_definitions(-DHAVE_ISO_FORTRAN_ENV)
endif()
CHECK_Fortran_SOURCE_COMPILES(
  "flush(5); end"
  HAVE_FLUSH_STMT
  SRC_EXT f90
  )
if(HAVE_FLUSH_STMT)
  add_definitions(-DHAVE_FLUSH_STMT)
endif()

#----------------------------------------------------------------------------
# Find MPI and set some flags so that FC and CC can point to gfortran and gcc
#----------------------------------------------------------------------------
# If the user passes FC=mpifort etc. check and prefer that location
get_filename_component( FTN_COMPILER_NAME "${CMAKE_Fortran_COMPILER}"
  NAME )
get_filename_component( C_COMPILER_NAME "${CMAKE_C_COMPILER}"
  NAME )
get_filename_component( FTN_COMPILER_DIR "${CMAKE_Fortran_COMPILER}"
  REALPATH )
get_filename_component( C_COMPILER_DIR "${CMAKE_C_COMPILER}"
  REALPATH )

if (FTN_COMPILER_NAME MATCHES "^[mM][pP][iI]")
  set (MPI_Fortran_COMPILER "${CMAKE_Fortran_COMPILER}")
endif()
if (C_COMPILER_NAME MATCHES "^[mM][pP][iI]")
  set (MPI_C_COMPILER "${CMAKE_C_COMPILER}")
endif()

find_package( MPI )

if ( (NOT MPI_C_FOUND) OR (NOT MPI_Fortran_FOUND) OR (NOT MPIEXEC))
  find_program (MY_MPI_EXEC NAMES mpirun mpiexec lamexec srun
    HINTS "${FTN_COMPILER_DIR}" "${C_COMPILER_DIR}"
    PATH_SUFFIXES bin)
  set ( MPI_HOME "${MPI_HOME}" "${MY_MPI_EXEC}" "${MY_MPI_EXEC}/.." )
  find_package( MPI REQUIRED )
endif()
list(REMOVE_DUPLICATES MPI_Fortran_INCLUDE_PATH)

# Test for consistent MPI environment
if (NOT MPIEXEC)
  message ( ERROR "CMake failed to find `mpiexec` or similar. If MPI is installed
on your system then please point CMake to the desired MPI runtime
using `cmake-gui` or `ccmake`. Otherwise, we will proceed with a
serial build. If you beleive you are getting this message in error,
please report this at https://github.com/sfilippone/psblas3/issues.")
endif()

get_filename_component(MPIEXEC_RELATIVE_LOC "${MPIEXEC}"
  PROGRAM)
get_filename_component(MPIEXEC_ABS_LOC "${MPIEXEC_RELATIVE_LOC}"
  REALPATH)
get_filename_component(MPIEXEC_DIR "${MPIEXEC_ABS_LOC}"
  DIRECTORY)

get_filename_component(MPICC_RELATIVE_LOC "${MPI_C_COMPILER}"
  PROGRAM)
get_filename_component(MPICC_ABS_LOC "${MPICC_RELATIVE_LOC}"
  REALPATH)
get_filename_component(MPICC_DIR "${MPICC_ABS_LOC}"
  DIRECTORY)

get_filename_component(MPIFC_RELATIVE_LOC "${MPI_Fortran_COMPILER}"
  PROGRAM)
get_filename_component(MPIFC_ABS_LOC "${MPIFC_RELATIVE_LOC}"
  REALPATH)
get_filename_component(MPIFC_DIR "${MPIFC_ABS_LOC}"
  DIRECTORY)

if ((MPIEXEC_DIR STREQUAL MPICC_DIR) AND (MPIEXEC_DIR STREQUAL MPIFC_DIR))
  message ( STATUS "MPI runtime and compile time environments appear to be consistent")
else()
  message ( WARNING "MPIEXEC is in \"${MPIEXEC_DIR},\"
which differs from the location of MPICC and/or MPIFC which are in
\"${MPICC_DIR}\" and \"${MPIFC_DIR},\" respectively.
This is likely indicative of an inconsistent MPI environment")
endif()

#-----------------------------------------------
# Work around an issue present on fedora systems
#-----------------------------------------------
if( (MPI_C_LINK_FLAGS MATCHES "noexecstack") OR (MPI_Fortran_LINK_FLAGS MATCHES "noexecstack") )
  message ( WARNING
"The `noexecstack` linker flag was found in the MPI_<lang>_LINK_FLAGS variable. This is
known to cause segmentation faults for some Fortran codes. See, e.g.,
https://gcc.gnu.org/bugzilla/show_bug.cgi?id=71729 or
https://github.com/sourceryinstitute/OpenCoarrays/issues/317.

`noexecstack` is being replaced with `execstack`"
    )
  string(REPLACE "noexecstack"
    "execstack" MPI_C_LINK_FLAGS_FIXED ${MPI_C_LINK_FLAGS})
  string(REPLACE "noexecstack"
    "execstack" MPI_Fortran_LINK_FLAGS_FIXED ${MPI_Fortran_LINK_FLAGS})
  set(MPI_C_LINK_FLAGS "${MPI_C_LINK_FLAGS_FIXED}" CACHE STRING
    "MPI C linking flags" FORCE)
  set(MPI_Fortran_LINK_FLAGS "${MPI_Fortran_LINK_FLAGS_FIXED}" CACHE STRING
    "MPI Fortran linking flags" FORCE)
endif()

#--------------------------------------------------------
# Make sure a simple "hello world" C mpi program compiles
#--------------------------------------------------------
set(OLD_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
set(CMAKE_REQUIRED_FLAGS ${MPI_C_COMPILE_FLAGS} ${MPI_C_LINK_FLAGS})
set(OLD_INCLUDES ${CMAKE_REQUIRED_INCLUDES})
set(CMAKE_REQUIRED_INCLUDES ${MPI_C_INCLUDE_PATH})
set(OLD_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES})
set(CMAKE_REQUIRED_LIBRARIES ${MPI_C_LIBRARIES})
include (CheckCSourceCompiles)
CHECK_C_SOURCE_COMPILES("
#include <mpi.h>
#include <stdio.h>
int main(int argc, char** argv) {
  MPI_Init(NULL, NULL);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);
  printf('Hello world from processor %s, rank %d out of %d processors',
         processor_name, world_rank, world_size);
  MPI_Finalize();
}"
MPI_C_COMPILES)
set(CMAKE_REQUIRED_FLAGS ${OLD_REQUIRED_FLAGS})
set(CMAKE_REQUIRED_INCLUDES ${OLD_INCLUDES})
set(CMAKE_REQUIRED_LIBRARIES ${OLD_LIBRARIES})
unset(OLD_REQUIRED_FLAGS)
unset(OLD_INCLUDES)
unset(OLD_LIBRARIES)

if (NOT MPI_C_COMPILES)
  message(FATAL_ERROR "MPI_C is missing! "
    "Try setting MPI_C_COMPILER to the appropriate C compiler wrapper script and reconfigure. "
    "i.e., `cmake -DMPI_C_COMPILER=/path/to/mpicc ..` or set it by editing the cache using "
    "cmake-gui or ccmake."
    )
endif()

#--------------------------------------------------------------
# Make sure a simple "hello world" Fortran mpi program compiles
# Try using mpi.mod first then fall back on includ 'mpif.h'
#--------------------------------------------------------------
set(OLD_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
set(CMAKE_REQUIRED_FLAGS "-ffree-form" ${MPI_Fortran_COMPILE_FLAGS} ${MPI_Fortran_LINK_FLAGS})
set(OLD_INCLUDES ${CMAKE_REQUIRED_INCLUDES})
set(CMAKE_REQUIRED_INCLUDES ${MPI_Fortran_INCLUDE_PATH})
set(OLD_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES})
set(CMAKE_REQUIRED_LIBRARIES ${MPI_Fortran_LIBRARIES})
include (CheckFortranSourceCompiles)
CHECK_Fortran_SOURCE_COMPILES("
program mpi_hello
use mpi
implicit none
integer :: ierr, mpi_world_size, mpi_world_rank, res_len
character*(MPI_MAX_PROCESSOR_NAME) :: proc
call mpi_init(ierr)
call mpi_comm_size(MPI_COMM_WORLD,mpi_world_size,ierr)
call mpi_comm_rank(MPI_COMM_WORLD,mpi_world_rank,ierr)
call mpi_get_processor_name(proc,res_len,ierr)
write(*,*) 'Hello from processor ', trim(proc), ' rank ', mpi_world_rank, ' out of ', mpi_world_size, '.'
call mpi_finalize(ierr)
end program
"
MPI_Fortran_MODULE_COMPILES)
set(CMAKE_REQUIRED_FLAGS ${OLD_REQUIRED_FLAGS})
set(CMAKE_REQUIRED_INCLUDES ${OLD_INCLUDES})
set(CMAKE_REQUIRED_LIBRARIES ${OLD_LIBRARIES})
unset(OLD_REQUIRED_FLAGS)
unset(OLD_INCLUDES)
unset(OLD_LIBRARIES)

#--------------------------------
# If that failed try using mpif.h
#--------------------------------
set(OLD_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
set(CMAKE_REQUIRED_FLAGS "-ffree-form" ${MPI_Fortran_COMPILE_FLAGS} ${MPI_Fortran_LINK_FLAGS})
set(OLD_INCLUDES ${CMAKE_REQUIRED_INCLUDES})
set(CMAKE_REQUIRED_INCLUDES ${MPI_Fortran_INCLUDE_PATH})
set(OLD_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES})
set(CMAKE_REQUIRED_LIBRARIES ${MPI_Fortran_LIBRARIES})
include (CheckFortranSourceCompiles)
CHECK_Fortran_SOURCE_COMPILES("
program mpi_hello
implicit none
include 'mpif.h'
integer :: ierr, mpi_world_size, mpi_world_rank, res_len
character*(MPI_MAX_PROCESSOR_NAME) :: proc
call mpi_init(ierr)
call mpi_comm_size(MPI_COMM_WORLD,mpi_world_size,ierr)
call mpi_comm_rank(MPI_COMM_WORLD,mpi_world_rank,ierr)
call mpi_get_processor_name(proc,res_len,ierr)
write(*,*) 'Hello from processor ', trim(proc), ' rank ', mpi_world_rank, ' out of ', mpi_world_size, '.'
call mpi_finalize(ierr)
end program
"
  MPI_Fortran_INCLUDE_COMPILES)
set(CMAKE_REQUIRED_FLAGS ${OLD_REQUIRED_FLAGS})
set(CMAKE_REQUIRED_INCLUDES ${OLD_INCLUDES})
set(CMAKE_REQUIRED_LIBRARIES ${OLD_LIBRARIES})
unset(OLD_REQUIRED_FLAGS)
unset(OLD_INCLUDES)
unset(OLD_LIBRARIES)

if ( (NOT MPI_Fortran_MODULE_COMPILES) AND (NOT MPI_Fortran_INCLUDE_COMPILES) )
  message ( WARNING "It appears that the Fortran MPI compiler is not working. ")
endif()

if ( MPI_Fortran_MODULE_COMPILES )
  add_definitions(-DMPI_WORKING_MODULE)
else()
  message ( WARNING "It appears that MPI was built with a different Fortran compiler. "
    "It is possible that this may cause unpredictable behavior. The build will continue "
    "using `mpif.h` BUT please report any suspicious behavior to the ${CMAKE_PROJECT_NAME} "
    "developers."
    )
endif()

#----------------
# Setup MPI flags
#----------------
set(CMAKE_C_COMPILE_FLAGS ${CMAKE_C_COMPILE_FLAGS} ${MPI_C_COMPILE_FLAGS})
set(CMAKE_C_LINK_FLAGS ${CMAKE_C_LINK_FLAGS} ${MPI_C_LINK_FLAGS})
set(CMAKE_Fortran_COMPILE_FLAGS ${CMAKE_Fortran_COMPILE_FLAGS} ${MPI_Fortran_COMPILE_FLAGS})
set(CMAKE_Fortran_LINK_FLAGS ${CMAKE_Fortran_LINK_FLAGS} ${MPI_Fortran_LINK_FLAGS})
include_directories(BEFORE ${MPI_C_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH})

#---------------------------------------------------
# Use standardized GNU install directory conventions
#---------------------------------------------------
include(GNUInstallDirs)
set(mod_dir_tail "${${CMAKE_PROJECT_NAME}_dist_string}_${CMAKE_Fortran_COMPILER_ID}-${CMAKE_Fortran_COMPILER_VERSION}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/${${CMAKE_PROJECT_NAME}_dist_string}-tests")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

enable_testing()


#------------------------------------------------------------------------------
# Add custom properties on targets for controling number of ranks during tests
#------------------------------------------------------------------------------
define_property(TARGET
  PROPERTY MIN_RANKS
  BRIEF_DOCS "Minimum allowable ranks for the test <integer>"
  FULL_DOCS "Property to mark executable targets run as tests that they require at least <MIN_RANKS> ranks to run"
  )

define_property(TARGET
  PROPERTY POWER_2_RANKS
  BRIEF_DOCS "True if test must be run with a power of 2 ranks (T/F)"
  FULL_DOCS "Property to mark executable targets run as tests that they require 2^n ranks."
  )

#-----------------------------------------------------
# Publicize installed location to other CMake projects
#-----------------------------------------------------
install(EXPORT ${CMAKE_PROJECT_NAME}-targets
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake"
)
include(CMakePackageConfigHelpers) # standard CMake module
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_PROJECT_NAME}ConfigVersion.cmake"
  VERSION "${psblas_VERSION}"
  COMPATIBILITY SameMajorVersion
  )

configure_file("${CMAKE_SOURCE_DIR}/cmake/pkg/${CMAKE_PROJECT_NAME}Config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${CMAKE_PROJECT_NAME}Config.cmake" @ONLY)

install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${CMAKE_PROJECT_NAME}Config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_PROJECT_NAME}ConfigVersion.cmake"
  DESTINATION
    "${CMAKE_INSTALL_LIBDIR}/cmake/psblas"
)


#------------------------------------------
# Add portable unistall command to makefile
#------------------------------------------
# Adapted from the CMake Wiki FAQ
configure_file ( "${CMAKE_SOURCE_DIR}/cmake/uninstall.cmake.in" "${CMAKE_BINARY_DIR}/uninstall.cmake"
    @ONLY)

add_custom_target ( uninstall
  COMMAND ${CMAKE_COMMAND} -P "${CMAKE_BINARY_DIR}/uninstall.cmake" )

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure)
# See JSON-Fortran's CMakeLists.txt file to find out how to get the check target to depend
# on the test executables

#---------------------------------------------------------------------------------------
# Define macro for adding CAF tests, and ensuring proper flags are passed to MPI runtime
#---------------------------------------------------------------------------------------

# Determine if we're using Open MPI
execute_process(COMMAND ${MPIEXEC} --version
  OUTPUT_VARIABLE mpi_version_out)
if (mpi_version_out MATCHES "[Oo]pen[ -][Mm][Pp][Ii]")
  message( STATUS "OpenMPI detected")
  set ( openmpi true )
endif ()

include( ProcessorCount )
ProcessorCount(N_CPU)
function(add_mpi_test name num_mpi_rank test_target)
  # Function to add MPI tests.
  if(TARGET ${test_target})
    get_target_property(min_test_ranks ${test_target} MIN_RANKS)
  elseif(TARGET build_${test_target})
    get_target_property(min_test_ranks build_${test_target} MIN_RANKS)
  endif()
  if(min_test_ranks)
    if(num_mpi_rank LESS min_test_ranks)
      message( FATAL_ERROR "Test ${name} requires ${min_test_ranks} but was only given ${num_mpi_ranks}" )
    endif()
  endif()
  if ( ((N_CPU LESS num_mpi_rank) OR (N_CPU EQUAL 0)) )
    message(STATUS "Test ${name} is oversubscribed: ${num_mpi_rank} CAF ranks requested with ${N_CPU} system processor available.")
    if ( openmpi )
      if (min_test_ranks)
	set( num_mpi_rank ${min_test_ranks} )
      elseif ( N_CPU LESS 2 )
        set( num_mpi_rank 2 )
      endif()
      set (test_parameters --oversubscribe)
      message( STATUS "Open-MPI back end detected, passing --oversubscribe for oversubscribed test, ${name}, with ${num_mpi_rank} ranks." )
    endif()
  endif()
  set(test_parameters -np ${num_mpi_rank} ${test_parameters})
  add_test(NAME ${name} COMMAND "bash" "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/mpirun" ${test_parameters} "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${test_target}")
  set_property(TEST ${name} PROPERTY PASS_REGULAR_EXPRESSION "Test passed.")
endfunction(add_mpi_test)

#-----------------
# Add PSBLAS tests
#-----------------

# Unit tests targeting each function, argument, and branch of code
# add_mpi_test(initialize_mpi 2 initialize_mpi)

add_definitions(
  -DHAVE_METIS -DHAVE_LAPACK -DHAVE_MOLD -DHAVE_EXTENDS_TYPE_OF -DHAVE_SAME_TYPE_AS -DHAVE_FINAL -DMPI_MOD
)
# From ./configure CDEFINES section
add_definitions(
  -DHAVE_METIS_
)

set(directory_list base prec krylov util test) # TODO: add cbind
foreach(directory ${directory_list})
  add_subdirectory("${directory}")
endforeach()
#add_subdirectory(cbind)
