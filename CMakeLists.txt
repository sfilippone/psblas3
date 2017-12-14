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
install(EXPORT PSBLAS_Targets
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

add_definitions(
  -DHAVE_METIS -DHAVE_LAPACK -DHAVE_MOLD -DHAVE_EXTENDS_TYPE_OF -DHAVE_SAME_TYPE_AS -DHAVE_FINAL -DHAVE_ISO_FORTRAN_ENV -DHAVE_FLUSH_STMT -DHAVE_VOLATILE -DMPI_MOD -DSERIAL_MPI
)
# From ./configure CDEFINES section
add_definitions(
  -DHAVE_METIS_  -DLowerUnderscore -DPtr64Bits
)

add_subdirectory(base)
