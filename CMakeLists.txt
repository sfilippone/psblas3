cmake_minimum_required(VERSION 3.2)

# Set the type/configuration of build to perform
set ( CMAKE_CONFIGURATION_TYPES "Debug" "Release" "MinSizeRel" "RelWithDebInfo" "CodeCoverage" )
set ( CMAKE_BUILD_TYPE "Release"
  CACHE STRING "Select which configuration to build." )
set_property ( CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${CMAKE_CONFIGURATION_TYPES} )

# Add option and check environment to determine if developer tests should be run
if($ENV{PSBLAS_DEVELOPER})
  option(MPI_RUN_DEVELOPER_TESTS "Run tests intended only for developers" ON)
else()
  option(MPI_RUN_DEVELOPER_TESTS "Run tests intended only for developers" OFF)
endif()
mark_as_advanced(MPI_RUN_DEVELOPER_TESTS)

if( NOT DEFINED ENV{PSBLAS_DEVELOPER})
  set ( ENV{PSBLAS_DEVELOPER} FALSE )
endif()

# Name project and specify source languages
# Parse version from .VERSION file so that more info can be added and easier to get from scripts
file(STRINGS ".VERSION" first_line
  LIMIT_COUNT 1
  )

string(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+"
  PSBLAS_Version "${first_line}")

if((NOT (PSBLAS_Version MATCHES "[0-9]+\\.[0-9]+\\.[0-9]+")) AND (EXISTS "${CMAKE_SOURCE_DIR}/.git"))
  message( STATUS "Build from git repository detected")
  find_package(Git)
  if(GIT_FOUND)
    execute_process(COMMAND "${GIT_EXECUTABLE}" describe --abbrev=0
      WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
      RESULT_VARIABLE git_status
      OUTPUT_VARIABLE git_output
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    if((git_status STREQUAL "0") AND (git_output MATCHES "[0-9]+\\.[0-9]+\\.[0-9]+"))
      set(PSBLAS_Version "${git_output}")
    endif()
    execute_process(COMMAND "${GIT_EXECUTABLE}" describe --always
      WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
      RESULT_VARIABLE git_status
      OUTPUT_VARIABLE full_git_describe
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT (git_status STREQUAL "0"))
      set(full_git_describe NOTFOUND)
    endif()
    # Create a source distribution target using git archive
    # e.g., `make dist` will package a release using current git state
    add_custom_target(dist # OUTPUT "${CMAKE_BINARY_DIR}/${_OC_stem_name}.tar.gz"
      COMMAND "${CMAKE_COMMAND}" -P "${CMAKE_SOURCE_DIR}/cmake/makeDist.cmake" "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}"
      COMMENT "Creating source release asset, ${_OC_stem_name}.tar.gz, from ${git_full_describe} using the `git archive` command."
      VERBATIM)
  else()
    message( WARNING "Could not find git executable!")
  endif()
endif()

if(NOT (PSBLAS_Version MATCHES "[0-9]+\\.[0-9]+\\.[0-9]+"))
  message( WARNING "Could not extract version from git, falling back on .VERSION, line 3.")
  file(STRINGS ".VERSION" PSBLAS_Version
  REGEX "[0-9]+\\.[0-9]+\\.[0-9]+"
  )
endif()

if(NOT full_git_describe)
  set(full_git_describe ${PSBLAS_Version})
endif()

# Strip leading "v" character from PSBLAS version tags
string(REPLACE "v" "" PSBLAS_Version "${PSBLAS_Version}")

project(psblas VERSION "${PSBLAS_Version}" LANGUAGES C Fortran)

message( STATUS "Building PSBLAS version: ${full_git_describe}" )
set(PSBLAS_dist_string "PSBLAS-${full_git_describe}")

#Print an error message on an attempt to build inside the source directory tree:
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
    "Please create a directory outside the PSBLAS source tree and build under that outside directory "
    "in a manner such as\n"
    "  $ mkdir build\n"
    "  $ cd build\n"
    "  $ CC=gcc FC=gfortran cmake -DBUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/path/to/install/dir /path/to/psblas3/src/dir \n"
    "\nsubstituting the appropriate syntax for your shell (the above line assumes the bash shell)."
    )
endif()

#Report untested Fortran compiler unless explicitly directed to build all examples.
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU" )
  set(gfortran_compiler true)
  set ( CMAKE_C_FLAGS_CODECOVERAGE "-fprofile-arcs -ftest-coverage -O0"
    CACHE STRING "Code coverage C compiler flags")
  set ( CMAKE_Fortran_FLAGS_CODECOVERAGE "-fprofile-arcs -ftest-coverage -O0"
    CACHE STRING "Code coverage Fortran compiler flags")
else()
  message(WARNING
    "\n"
    "Attempting to build with untested Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}. "
    "Please report any failures at https://github.com/sfilippone/psblas3\n\n"
  )
endif()

#-----------------------------------------------------------------
# Set CMAKE_Fortran_COMPILER_VERSION if CMake doesn't do it for us
#-----------------------------------------------------------------
if (  NOT CMAKE_Fortran_COMPILER_VERSION )
  if ( NOT (CMAKE_VERSION VERSION_LESS 3.3.1) )
    message( AUTHOR_WARNING
     "CMake ${CMAKE_VERSION} should know about Fortran compiler versions but is missing CMAKE_Fortran_COMPILER_VERSION variable."
    )
  endif()
  # No CMAKE_Fortran_COMPILER_VERSION set, build our own
  # Try extracting it directly from ISO_FORTRAN_ENV's compiler_version
  # Write program for introspection
  file( WRITE "${CMAKE_BINARY_DIR}/get_compiler_ver.f90"
    "program main
     use iso_fortran_env, only: compiler_version, output_unit
     write(output_unit,'(a)') compiler_version()
end program"
  )
  try_run( PROG_RAN COMPILE_SUCCESS
    "${CMAKE_BINARY_DIR}" "${CMAKE_BINARY_DIR}/get_compiler_ver.f90"
    RUN_OUTPUT_VARIABLE VER_STRING
  )
  if ( COMPILE_SUCCESS )
    string( REGEX MATCH "[0-9]+\\.[0-9]+(\\.[0-9]+)?"
      DETECTED_VER "${VER_STRING}"
    )
    message( STATUS "Detected Fortran compiler as ${VER_STRING}" )
    message( STATUS "Extracted version number: ${DETECTED_VER}" )
  endif()
  if( ( NOT COMPILE_SUCCESS ) OR ( NOT DETECTED_VER ) )
    message( WARNING "Could not reliably detect Fortran compiler version. We'll infer it from
the C compiler if it matches the Fortran compiler ID." )
  endif()
  if( "${CMAKE_C_COMPILER_ID}" MATCHES "${CMAKE_Fortran_COMPILER_ID}" )
    set( DETECTED_VER "${CMAKE_C_COMPILER_VERSION}" )
  else()
    message( FATAL_ERROR "Exhausted all possible means of detecting the Fortran compiler version, cannot proceed!" )
  endif()
  set( CMAKE_Fortran_COMPILER_VERSION "${DETECTED_VER}" )
endif()

if(gfortran_compiler)
  set(OLD_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
  set(CMAKE_REQUIRED_FLAGS "-fcoarray=single -ffree-form")
endif()
include(CheckFortranSourceCompiles)
CHECK_Fortran_SOURCE_COMPILES("
  program main
    implicit none
    integer :: i
    i = this_image()
  end program
" Check_Simple_Coarray_Fortran_Source_Compiles)
if(gfortran_compiler)
  set (CMAKE_REQUIRED_FLAGS ${OLD_REQUIRED_FLAGS})
  unset(OLD_REQUIRED_FLAGS)
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
  # Get default install location of MPICH from install.sh
  execute_process( COMMAND "./install.sh" -P mpich
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    OUTPUT_VARIABLE DEFAULT_MPICH_INSTALL_LOC
    OUTPUT_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  find_program (MY_MPI_EXEC NAMES mpirun mpiexec lamexec srun
    PATHS "${DEFAULT_MPICH_INSTALL_LOC}" ENV PATH
    HINTS "${FTN_COMPILER_DIR}" "${C_COMPILER_DIR}"
    PATH_SUFFIXES bin)
  set ( MPI_HOME "${MPI_HOME}" "${MY_MPI_EXEC}" "${MY_MPI_EXEC}/.." )
  find_package( MPI REQUIRED )
endif()
list(REMOVE_DUPLICATES MPI_Fortran_INCLUDE_PATH)

# Test for consistent MPI environment
if (NOT MPIEXEC)
  message ( ERROR "CMake failed to find `mpiexec` or similar. If building with `./install.sh` please
report this bug to the PSBLAS developers at
https://github.com/sfilippone/psblas3/issues, otherwise point CMake
to the desired MPI runtime.")
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
This is likely indicative of a problem. If building with `./install.sh` please report
this to the PSBLAS developers by filing a new issue at:
https://github.com/sourceryinstitute/PSBLAS/issues/new")
endif()

#-----------------------------------------------
# Work around bug #317 present on fedora systems
#-----------------------------------------------
if( (MPI_C_LINK_FLAGS MATCHES "noexecstack") OR (MPI_Fortran_LINK_FLAGS MATCHES "noexecstack") )
  message ( WARNING
"The `noexecstack` linker flag was found in the MPI_<lang>_LINK_FLAGS variable. This is
known to cause segmentation faults for some Fortran codes. See, e.g.,
https://gcc.gnu.org/bugzilla/show_bug.cgi?id=71729 or
https://github.com/sourceryinstitute/PSBLAS/issues/317.

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
    "using `mpif.h` BUT please report any suspicious behavior to the PSBLAS "
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
set(mod_dir_tail "${PSBLAS_dist_string}_${CMAKE_Fortran_COMPILER_ID}-${CMAKE_Fortran_COMPILER_VERSION}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/${PSBLAS_dist_string}-tests")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

#-----------------
# Install manpages
#-----------------
install(FILES "${CMAKE_SOURCE_DIR}/doc/man/man1/caf.1" "${CMAKE_SOURCE_DIR}/doc/man/man1/cafrun.1"
  DESTINATION "${CMAKE_INSTALL_MANDIR}/man1"
  PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ
  COMPONENT documentation)

#---------------------------------------------------
# Define macro for compiling with wrapper script
#---------------------------------------------------
function(caf_compile_executable target main_depend)
  foreach(includedir IN LISTS MPI_Fortran_INCLUDE_PATH)
    set(includes ${includes} -I ${includedir})
  endforeach()
  string(TOUPPER "${CMAKE_BUILD_TYPE}" build_type)
  separate_arguments(config_Fortran_flags UNIX_COMMAND "${CMAKE_Fortran_FLAGS_${build_type}}")
  get_directory_property( DirDefs DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" COMPILE_DEFINITIONS )
  set(localDefs "")
  foreach(d ${DirDefs})
    list(APPEND localDefs "-D${d}")
  endforeach()
  add_custom_command(OUTPUT "${target}"
    COMMAND "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/caf" ${includes} ${localDefs} ${config_Fortran_flags} -o "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${target}" "${CMAKE_CURRENT_SOURCE_DIR}/${main_depend}" ${ARGN}
    DEPENDS "${main_depend}" ${ARGN} caf_mpi_static
    VERBATIM
    )
  add_custom_target("build_${target}" ALL
    DEPENDS "${target}")
endfunction(caf_compile_executable)

enable_testing()

#--------------------------------------------------------
# Setup shellcheck if present for testing/linting scripts
#--------------------------------------------------------
find_program(SHELLCHECK_EXE shellcheck
  DOC "Path to shellcheck executable for linting scripts"
  )
if (MPI_RUN_DEVELOPER_TESTS OR $ENV{PSBLAS_DEVELOPER})
  if(NOT SHELLCHECK_EXE)
    message( AUTHOR_WARNING "PSBLAS developers should install shellcheck to test/lint all shell scripts.
    See https://github.com/koalaman/shellcheck#installing for info on obtaining shellcheck.")
  endif()
endif()

function(lint_script script_dir script_name)
  if (SHELLCHECK_EXE)
    add_test(NAME "shellcheck:${script_name}"
      COMMAND ${SHELLCHECK_EXE} -x "${script_dir}/${script_name}"
      WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
  elseif (MPI_RUN_DEVELOPER_TESTS OR $ENV{PSBLAS_DEVELOPER})
    message( AUTHOR_WARNING "test: shellcheck:${script_name} not run because shellcheck not installed." )
  endif()
endfunction()

#-----------------------------------------------
# Setup script style testing & enforcement macro
#-----------------------------------------------

find_program(style_pl style.pl "${CMAKE_SOURCE_DIR}/developer-scripts/")
function(check_script_style script_full_path)
  if(style_pl)
    add_test(NAME "style:${script_full_path}"
      COMMAND "${style_pl}" "${script_full_path}")
  endif()
endfunction()

#------------------------------------------------------------------------------
# Add custom properties on targets for controling number of images during tests
#------------------------------------------------------------------------------
define_property(TARGET
  PROPERTY MIN_IMAGES
  BRIEF_DOCS "Minimum allowable images for the test <integer>"
  FULL_DOCS "Property to mark executable targets run as tests that they require at least <MIN_IMAGES> images to run"
  )

define_property(TARGET
  PROPERTY POWER_2_IMGS
  BRIEF_DOCS "True if test must be run with a power of 2 images (T/F)"
  FULL_DOCS "Property to mark executable targets run as tests that they require 2^n images."
  )


#-------------------------------
# Recurse into the src directory
#-------------------------------
#include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/src)

#add_subdirectory(src)

#-----------------------------------------------------
# Publicize installed location to other CMake projects
#-----------------------------------------------------
install(EXPORT PSBLASTargets
  NAMESPACE
    PSBLAS::
  DESTINATION
    "${CMAKE_INSTALL_LIBDIR}/cmake/psblas"
)
include(CMakePackageConfigHelpers) # standard CMake module
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/PSBLASConfigVersion.cmake"
  VERSION "${psblas_VERSION}"
  COMPATIBILITY AnyNewerVersion
)
configure_file("${CMAKE_SOURCE_DIR}/cmake/pkg/PSBLASConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/PSBLASConfig.cmake" @ONLY)

install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/PSBLASConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/PSBLASConfigVersion.cmake"
  DESTINATION
    "${CMAKE_INSTALL_LIBDIR}/cmake/psblas"
)

add_library(PSBLAS INTERFACE)
#target_compile_options(PSBLAS INTERFACE -fcoarray=lib)
#target_link_libraries(PSBLAS INTERFACE caf_mpi)

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
    get_target_property(min_test_ranks ${test_target} MIN_IMAGES)
  elseif(TARGET build_${test_target})
    get_target_property(min_test_ranks build_${test_target} MIN_IMAGES)
  endif()
  if(min_test_ranks)
    if(num_mpi_rank LESS min_test_ranks)
      message( FATAL_ERROR "Test ${name} requires ${min_test_ranks} but was only given ${num_mpi_images}" )
    endif()
  endif()
  if ( ((N_CPU LESS num_mpi_rank) OR (N_CPU EQUAL 0)) )
    message(STATUS "Test ${name} is oversubscribed: ${num_mpi_rank} CAF images requested with ${N_CPU} system processor available.")
    if ( openmpi )
      if (min_test_ranks)
	set( num_mpi_rank ${min_test_ranks} )
      elseif ( N_CPU LESS 2 )
        set( num_mpi_rank 2 )
      endif()
      set (test_parameters --oversubscribe)
      message( STATUS "Open-MPI back end detected, passing --oversubscribe for oversubscribed test, ${name}, with ${num_mpi_rank} ranks/images." )
    endif()
  endif()
  set(test_parameters -np ${num_mpi_rank} ${test_parameters})
  add_test(NAME ${name} COMMAND "bash" "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/mpirun" ${test_parameters} "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${test_target}")
  set_property(TEST ${name} PROPERTY PASS_REGULAR_EXPRESSION "Test passed.")
endfunction(add_mpi_test)

#-----------------
# Add PSBLAS tests
#-----------------

if (MPI_RUN_DEVELOPER_TESTS OR $ENV{PSBLAS_DEVELOPER})
  message ( STATUS "Running Developer tests is enabled." )
endif()

# Unit tests targeting each function, argument, and branch of code
# add_mpi_test(initialize_mpi 2 initialize_mpi)

#include(cmake/AddInstallationScriptTest.cmake )
#add_installation_script_test(installation-scripts.sh src/tests/installation/)

add_library(base
base/internals/psi_desc_impl.f90
base/internals/psi_crea_bnd_elem.f90
base/internals/psi_bld_tmphalo.f90
base/internals/psi_list_search.f90
base/internals/psi_crea_ovr_elem.f90
base/internals/psi_compute_size.f90
base/internals/psi_dl_check.f90
base/internals/psi_exist_ovr_elem.f90
base/internals/psi_srtlist.f90
base/internals/psi_bld_tmpovrl.f90
base/internals/psi_crea_index.f90
base/internals/psi_sort_dl.f90
base/tools/psb_iallc.f90
base/tools/psb_c_map.f90
base/tools/psb_dspalloc.f90
base/tools/psb_cspasb.f90
base/tools/psb_dins.f90
base/tools/psb_iasb.f90
base/tools/psb_sspfree.f90
base/tools/psb_cd_lstext.f90
base/tools/psb_dsprn.f90
base/tools/psb_cins.f90
base/tools/psb_csprn.f90
base/tools/psb_ifree.f90
base/tools/psb_sspasb.f90
base/tools/psb_z_map.f90
base/tools/psb_zspins.f90
base/tools/psb_sins.f90
base/tools/psb_cdins.f90
base/tools/psb_cspfree.f90
base/tools/psb_cdalv.f90
base/tools/psb_zallc.f90
base/tools/psb_cspins.f90
base/tools/psb_dallc.f90
base/tools/psb_cdrep.f90
base/tools/psb_dspfree.f90
base/tools/psb_zspfree.f90
base/tools/psb_get_overlap.f90
base/tools/psb_cdals.f90
base/tools/psb_dfree.f90
base/tools/psb_dasb.f90
base/tools/psb_zfree.f90
base/tools/psb_zins.f90
base/tools/psb_cd_switch_ovl_indxmap.f90
base/tools/psb_sasb.f90
base/tools/psb_callc.f90
base/tools/psb_cd_reinit.f90
base/tools/psb_cd_set_bld.f90
base/tools/psb_s_map.f90
base/tools/psb_cdprt.f90
base/tools/psb_glob_to_loc.f90
base/tools/psb_zasb.f90
base/tools/psb_zspalloc.f90
base/tools/psb_cspalloc.f90
base/tools/psb_iins.f90
base/tools/psb_d_map.f90
base/tools/psb_cdren.f90
base/tools/psb_casb.f90
base/tools/psb_zsprn.f90
base/tools/psb_sspins.f90
base/tools/psb_zspasb.f90
base/tools/psb_loc_to_glob.f90
base/tools/psb_dspasb.f90
base/tools/psb_cdall.f90
base/tools/psb_sallc.f90
base/tools/psb_ssprn.f90
base/tools/psb_dspins.f90
base/tools/psb_sspalloc.f90
base/tools/psb_cd_inloc.f90
base/tools/psb_sfree.f90
base/tools/psb_cfree.f90
base/modules/psi_c_mod.f90
base/modules/psi_mod.f90
base/modules/aux/psb_string_mod.f90
base/modules/aux/psb_sort_mod.f90
base/modules/aux/psb_s_sort_mod.f90
base/modules/aux/psi_d_serial_mod.f90
base/modules/aux/psi_serial_mod.f90
base/modules/aux/psb_i_sort_mod.f90
base/modules/aux/psi_z_serial_mod.f90
base/modules/aux/psi_i_serial_mod.f90
base/modules/aux/psb_ip_reord_mod.f90
base/modules/aux/psb_c_sort_mod.f90
base/modules/aux/psb_hash_mod.f90
base/modules/aux/psb_d_sort_mod.f90
base/modules/aux/psi_s_serial_mod.f90
base/modules/aux/psi_c_serial_mod.f90
base/modules/aux/psb_z_sort_mod.f90
base/modules/psb_base_mod.f90
base/modules/psi_z_mod.f90
base/modules/desc/psb_indx_map_mod.f90
base/modules/desc/psb_repl_map_mod.f90
base/modules/desc/psb_desc_const_mod.f90
base/modules/desc/psb_gen_block_map_mod.f90
base/modules/desc/psb_hash_map_mod.f90
base/modules/desc/psb_glist_map_mod.f90
base/modules/desc/psb_list_map_mod.f90
base/modules/psi_d_mod.f90
base/modules/tools/psb_i_tools_mod.f90
base/modules/tools/psb_d_tools_mod.f90
base/modules/tools/psb_s_tools_mod.f90
base/modules/tools/psb_tools_mod.f90
base/modules/tools/psb_cd_tools_mod.f90
base/modules/tools/psb_z_tools_mod.f90
base/modules/tools/psb_c_tools_mod.f90
base/modules/psi_s_mod.f90
base/modules/comm/psb_d_linmap_mod.f90
base/modules/comm/psb_s_linmap_mod.f90
base/modules/comm/psb_base_linmap_mod.f90
base/modules/comm/psb_i_comm_mod.f90
base/modules/comm/psb_c_linmap_mod.f90
base/modules/comm/psb_c_comm_mod.f90
base/modules/comm/psb_z_comm_mod.f90
base/modules/comm/psb_z_linmap_mod.f90
base/modules/comm/psb_s_comm_mod.f90
base/modules/comm/psb_d_comm_mod.f90
base/modules/comm/psb_comm_mod.f90
base/modules/comm/psb_linmap_mod.f90
base/modules/psblas/psb_psblas_mod.f90
base/modules/psi_i_mod.f90
base/modules/serial/psb_s_csr_mat_mod.f90
base/modules/serial/psb_z_mat_mod.f90
base/modules/serial/psb_s_base_mat_mod.f90
base/modules/serial/psb_s_base_vect_mod.f90
base/modules/serial/psb_z_csc_mat_mod.f90
base/modules/serial/psb_d_csr_mat_mod.f90
base/modules/serial/psb_s_csc_mat_mod.f90
base/modules/serial/psb_i_base_vect_mod.f90
base/modules/serial/psb_z_base_vect_mod.f90
base/modules/serial/psb_c_csr_mat_mod.f90
base/modules/serial/psb_c_serial_mod.f90
base/modules/serial/psb_c_base_vect_mod.f90
base/modules/serial/psb_s_serial_mod.f90
base/modules/serial/psb_d_csc_mat_mod.f90
base/modules/serial/psb_vect_mod.f90
base/modules/serial/psb_base_mat_mod.f90
base/modules/serial/psb_c_csc_mat_mod.f90
base/modules/serial/psb_z_base_mat_mod.f90
base/modules/serial/psb_z_csr_mat_mod.f90
base/modules/serial/psb_c_mat_mod.f90
base/modules/serial/psb_d_base_mat_mod.f90
base/modules/serial/psb_d_serial_mod.f90
base/modules/serial/psb_c_base_mat_mod.f90
base/modules/serial/psb_d_base_vect_mod.f90
base/modules/serial/psb_serial_mod.f90
base/modules/serial/psb_s_mat_mod.f90
base/modules/serial/psb_z_serial_mod.f90
base/modules/serial/psb_d_mat_mod.f90
base/modules/serial/psb_mat_mod.f90
base/modules/error.f90
base/modules/psb_check_mod.f90
base/comm/psb_shalo.f90
base/comm/psb_dhalo.f90
base/comm/psb_igather.f90
base/comm/internals/psi_iovrl_restr.f90
base/comm/internals/psi_covrl_upd.f90
base/comm/internals/psi_covrl_save.f90
base/comm/internals/psi_sovrl_upd.f90
base/comm/internals/psi_dovrl_restr.f90
base/comm/internals/psi_zovrl_restr.f90
base/comm/internals/psi_dovrl_save.f90
base/comm/internals/psi_zovrl_save.f90
base/comm/internals/psi_sovrl_restr.f90
base/comm/internals/psi_zovrl_upd.f90
base/comm/internals/psi_sovrl_save.f90
base/comm/internals/psi_iovrl_upd.f90
base/comm/internals/psi_covrl_restr.f90
base/comm/internals/psi_iovrl_save.f90
base/comm/internals/psi_dovrl_upd.f90
base/comm/psb_iovrl.f90
base/comm/psb_cgather.f90
base/comm/psb_covrl.f90
base/comm/psb_zovrl.f90
base/comm/psb_ihalo.f90
base/comm/psb_dgather.f90
base/comm/psb_sovrl.f90
base/comm/psb_dovrl.f90
base/comm/psb_zgather.f90
base/comm/psb_zhalo.f90
base/comm/psb_sgather.f90
base/comm/psb_chalo.f90
base/psblas/psb_snrmi.f90
base/psblas/psb_dspmm.f90
base/psblas/psb_samax.f90
base/psblas/psb_sxdot.f90
base/psblas/psb_dspnrm1.f90
base/psblas/psb_sspnrm1.f90
base/psblas/psb_sspmm.f90
base/psblas/psb_dnrmi.f90
base/psblas/psb_zamax.f90
base/psblas/psb_cspmm.f90
base/psblas/psb_dspsm.f90
base/psblas/psb_casum.f90
base/psblas/psb_cnrmi.f90
base/psblas/psb_znrm2.f90
base/psblas/psb_zspsm.f90
base/psblas/psb_camax.f90
base/psblas/psb_dnrm2.f90
base/psblas/psb_znrmi.f90
base/psblas/psb_cspsm.f90
base/psblas/psb_sdot.f90
base/psblas/psb_dasum.f90
base/psblas/psb_zasum.f90
base/psblas/psb_damax.f90
base/psblas/psb_zspmm.f90
base/psblas/psb_cnrm2.f90
base/psblas/psb_ddot.f90
base/psblas/psb_caxpby.f90
base/psblas/psb_cdot.f90
base/psblas/psb_snrm2.f90
base/psblas/psb_sspsm.f90
base/psblas/psb_daxpby.f90
base/psblas/psb_zdot.f90
base/psblas/psb_zaxpby.f90
base/psblas/psb_zspnrm1.f90
base/psblas/psb_sasum.f90
base/psblas/psb_saxpby.f90
base/psblas/psb_cspnrm1.f90
base/serial/psb_zspspmm.f90
base/serial/psi_z_serial_impl.f90
base/serial/psb_dgelp.f90
base/serial/psb_zrwextd.f90
base/serial/sort/psb_c_hsort_impl.f90
base/serial/sort/psb_d_msort_impl.f90
base/serial/sort/psi_acx_mod.f90
base/serial/sort/psb_s_hsort_impl.f90
base/serial/sort/psb_i_isort_impl.f90
base/serial/sort/psb_z_qsort_impl.f90
base/serial/sort/psi_lcx_mod.f90
base/serial/sort/psb_s_isort_impl.f90
base/serial/sort/psb_i_msort_impl.f90
base/serial/sort/psi_alcx_mod.f90
base/serial/sort/psb_c_qsort_impl.f90
base/serial/sort/psb_c_msort_impl.f90
base/serial/sort/psb_z_hsort_impl.f90
base/serial/sort/psb_z_isort_impl.f90
base/serial/sort/psb_d_qsort_impl.f90
base/serial/sort/psb_s_qsort_impl.f90
base/serial/sort/psb_d_hsort_impl.f90
base/serial/sort/psb_c_isort_impl.f90
base/serial/sort/psb_i_hsort_impl.f90
base/serial/sort/psb_s_msort_impl.f90
base/serial/sort/psb_z_msort_impl.f90
base/serial/sort/psb_d_isort_impl.f90
base/serial/sort/psb_i_qsort_impl.f90
base/serial/psb_dgeprt.f90
base/serial/psb_sgelp.f90
base/serial/psb_zgeprt.f90
base/serial/psb_casum_s.f90
base/serial/psb_lsame.f90
base/serial/psb_dasum_s.f90
base/serial/impl/psb_z_coo_impl.f90
base/serial/impl/psb_z_csr_impl.f90
base/serial/impl/psb_d_csc_impl.f90
base/serial/impl/psb_base_mat_impl.f90
base/serial/impl/psb_d_coo_impl.f90
base/serial/impl/psb_z_csc_impl.f90
base/serial/impl/psb_d_csr_impl.f90
base/serial/impl/psb_s_coo_impl.f90
base/serial/impl/psb_c_csr_impl.f90
base/serial/impl/psb_s_csr_impl.f90
base/serial/impl/psb_s_csc_impl.f90
base/serial/impl/psb_c_csc_impl.f90
base/serial/impl/psb_c_coo_impl.f90
base/serial/smmp.f90
base/serial/psb_srwextd.f90
base/serial/psb_camax_s.f90
base/serial/psb_ssymbmm.f90
base/serial/psb_cgeprt.f90
base/serial/psb_zasum_s.f90
base/serial/psi_d_serial_impl.f90
base/serial/psb_znumbmm.f90
base/serial/psb_cspspmm.f90
base/serial/psb_sgeprt.f90
base/serial/psb_samax_s.f90
base/serial/psb_sasum_s.f90
base/serial/psi_i_serial_impl.f90
base/serial/psb_aspxpby.f90
base/serial/psb_zamax_s.f90
base/serial/psi_c_serial_impl.f90
base/serial/psb_sspspmm.f90
base/serial/psb_cnumbmm.f90
base/serial/psi_s_serial_impl.f90
base/serial/psb_dspspmm.f90
base/serial/psb_zsymbmm.f90
base/serial/psb_cgelp.f90
base/serial/psb_dnumbmm.f90
base/serial/psb_csymbmm.f90
base/serial/psb_zgelp.f90
base/serial/psb_dsymbmm.f90
base/serial/psb_damax_s.f90
base/serial/psb_snumbmm.f90
base/serial/psb_drwextd.f90
base/serial/psb_spge_dot.f90
base/serial/psb_spdot_srtd.f90
base/serial/psb_crwextd.f90
base/internals/psb_indx_map_fnd_owner.F90
base/internals/psi_desc_index.F90
base/internals/psi_fnd_owner.F90
base/internals/psi_extrct_dl.F90
base/tools/psb_icdasb.F90
base/tools/psb_dsphalo.F90
base/tools/psb_zcdbldext.F90
base/tools/psb_scdbldext.F90
base/tools/psb_zsphalo.F90
base/tools/psb_ccdbldext.F90
base/tools/psb_cdcpy.F90
base/tools/psb_ssphalo.F90
base/tools/psb_csphalo.F90
base/tools/psb_dcdbldext.F90
base/modules/psb_penv_mod.F90
base/modules/psi_comm_buffers_mod.F90
base/modules/psi_p2p_mod.F90
base/modules/desc/psb_desc_mod.F90
base/modules/psi_penv_mod.F90
base/modules/psb_error_impl.F90
base/modules/psb_realloc_mod.F90
base/modules/psb_error_mod.F90
base/modules/psblas/psb_s_psblas_mod.F90
base/modules/psblas/psb_d_psblas_mod.F90
base/modules/psblas/psb_c_psblas_mod.F90
base/modules/psblas/psb_z_psblas_mod.F90
base/modules/psi_reduce_mod.F90
base/modules/psb_const_mod.F90
base/modules/serial/psb_i_vect_mod.F90
base/modules/serial/psb_c_vect_mod.F90
base/modules/serial/psb_d_vect_mod.F90
base/modules/serial/psb_z_vect_mod.F90
base/modules/serial/psb_s_vect_mod.F90
base/modules/psi_bcast_mod.F90
base/comm/psb_zspgather.F90
base/comm/psb_dscatter.F90
base/comm/psb_zscatter.F90
base/comm/internals/psi_dswapdata.F90
base/comm/internals/psi_sswaptran.F90
base/comm/internals/psi_sswapdata.F90
base/comm/internals/psi_cswaptran.F90
base/comm/internals/psi_iswapdata.F90
base/comm/internals/psi_zswaptran.F90
base/comm/internals/psi_iswaptran.F90
base/comm/internals/psi_dswaptran.F90
base/comm/internals/psi_zswapdata.F90
base/comm/internals/psi_cswapdata.F90
base/comm/psb_cspgather.F90
base/comm/psb_iscatter.F90
base/comm/psb_sscatter.F90
base/comm/psb_sspgather.F90
base/comm/psb_cscatter.F90
base/comm/psb_dspgather.F90
base/serial/impl/psb_c_base_mat_impl.F90
base/serial/impl/psb_s_mat_impl.F90
base/serial/impl/psb_z_mat_impl.F90
base/serial/impl/psb_d_base_mat_impl.F90
base/serial/impl/psb_c_mat_impl.F90
base/serial/impl/psb_s_base_mat_impl.F90
base/serial/impl/psb_d_mat_impl.F90
base/serial/impl/psb_z_base_mat_impl.F90
base/modules/fakempi.c
base/modules/cutil.c
)
