# FDEFINES : -DHAVE_LAPACK -DHAVE_FINAL -DHAVE_ISO_FORTRAN_ENV -DHAVE_FLUSH_STMT -DHAVE_VOLATILE -DSERIAL_MPI -DMPI_MOD
# CDEFINES : -DLowerUnderscore -DPtr64Bits

#-----------------------------------
# Set oldest allowable CMake version
#-----------------------------------
cmake_minimum_required(VERSION 3.20)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

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
include(setVersion)
set_version(
  VERSION_VARIABLE PSBLAS_Version
  GIT_DESCRIBE_VAR full_git_describe
  CUSTOM_VERSION_FILE "${CMAKE_SOURCE_DIR}/.VERSION")
message( STATUS "Building PSBLAS1 version: ${full_git_describe}" )
#------------------------------------------
# Name project and specify source languages
#------------------------------------------
project(psblas
  VERSION "${PSBLAS_Version}"
  LANGUAGES C Fortran)

#--------------------------------------------------
# Set option to allow building against OpenCoarrays
#--------------------------------------------------
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  option(PSBLAS_USE_OpenCoarrays "Build enabling linkage to programs using OpenCoarrays" OFF)
endif()

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
include(CheckOutOfSourceBuild)

#----------------------------------------------------
# Define coverage flags and report untested compilers
#----------------------------------------------------
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU" )
  set(CMAKE_Fortran_COMPILER mpifort CACHE STRING "MPI Fortran compiler" FORCE)

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

#----------------------------------------------------
# Define -frecursive for GNU Fortran Compiler
#----------------------------------------------------
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU" )
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -frecursive")
    message(STATUS "GNU Fortran COMPILER ${CMAKE_Fortran_FLAGS};")
endif()


  message(STATUS "cmake flags? ${CMAKE_Fortran_FLAGS};")

#------------------------------------
# Fortran name mangling introspection
#------------------------------------
include("${CMAKE_CURRENT_LIST_DIR}/cmake/CapitalizeString.cmake")
#include(FortranCInterface)
#CapitalizeString(${FortranCInterface_GLOBAL__CASE} fc_case)
#message(STATUS "Name mangling capitalization: ${fc_case}")
#message(STATUS "Name mangling fortran global suffix underscore: ${FortranCInterface_GLOBAL__SUFFIX}")
#if(FortranCInterface_GLOBAL__SUFFIX STREQUAL "")
#  add_compile_options("-D${fc_case}Case")
#elseif(FortranCInterface_GLOBAL__SUFFIX STREQUAL "_")
#  add_compile_options("-D${fc_case}Underscore")
#elseif(FortranCInterface_GLOBAL__SUFFIX STREQUAL "__")
#  add_compile_options("-D${fc_case}DoubleUnderscore")
#else()
#  message( FATAL_ERROR "Fortran name mangling suffix, \'${FortranCInterface_GLOBAL__SUFFIX}\', unknown to PSBLAS")
#endif()


#  message(STATUS "win? ${WIN32};")
#if(TRUE)#NOT ${WIN32})
#previous check did not work if WIN32 is empty string
  #----------------------------------------------
  # Determine system endian-ness and pointer size
  #----------------------------------------------
#  include(TestBigEndian)
#  TEST_BIG_ENDIAN(IS_BIG_ENDIAN)
#  if(IS_BIG_ENDIAN)
#    message( STATUS "System appears to be big endian.")
#  else()
#    message( STATUS "System appears to be little endian.")
#    add_compile_options(-DLittleEndian)
#  endif()
#  include(CheckTypeSize)
#  CHECK_TYPE_SIZE("void *" VOID_P_SIZE LANGUAGE C)
#  if(${VOID_P_SIZE} EQUAL 8)
#    add_compile_options(-DPtr64Bits)
#  endif()
#  message(STATUS "Have 64bit pointers")
 
  
#endif()

message(STATUS "Using compiler ${CMAKE_C_COMPILER};")


# Set default values for IPK_SIZE and LPK_SIZE
set(DEFAULT_IPK_SIZE 4)
set(DEFAULT_LPK_SIZE 8)

# Allow user to override with command line definitions
if(NOT DEFINED CMAKE_PSB_IPK)
    set(CMAKE_PSB_IPK ${DEFAULT_IPK_SIZE} CACHE STRING "Size of IPK (default: 4)")
endif()

if(NOT DEFINED CMAKE_PSB_LPK)  
    set(CMAKE_PSB_LPK ${DEFAULT_LPK_SIZE} CACHE STRING "Size of LPK (default: 8)")
endif()

# Use the passed values
set(IPK_SIZE ${CMAKE_PSB_IPK})
set(LPK_SIZE ${CMAKE_PSB_LPK})			
# Define IPKDEF and LPKDEF based on the sizes
set(PSB_IPKDEF "#define PSB_IPK${IPK_SIZE}")
set(PSB_LPKDEF "#define PSB_LPK${LPK_SIZE}")

# Output the definitions for verification
message(STATUS "Using IPKDEF: ${PSB_IPKDEF}")
message(STATUS "Using LPKDEF: ${PSB_LPKDEF}")

#add_compile_options(-DPSB_IPK${IPK_SIZE})
#add_compile_options(-DPSB_LPK${LPK_SIZE})
# Add PSB_IPK/LPK flag only for fortran files.
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DPSB_IPK${IPK_SIZE}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DPSB_LPK${LPK_SIZE}")


  

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
  #add_compile_options(-DHAVE_MOVE_ALLOC)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DHAVE_MOVE_ALLOC")
  message(STATUS "-DHAVE_MOVE_ALLOC")
endif()
CHECK_Fortran_SOURCE_COMPILES(
  "integer, volatile :: i ; end"
  HAVE_VOLATILE
  SRC_EXT f90
  )
if(HAVE_VOLATILE)
  #add_compile_options(-DPSB_HAVE_VOLATILE)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DPSB_HAVE_VOLATILE")
  message(STATUS "-DPSB_HAVE_VOLATILE")
endif()
CHECK_Fortran_SOURCE_COMPILES(
  "use ISO_FORTRAN_ENV ; end"
  HAVE_ISO_FORTRAN_ENV
  SRC_EXT f90
  )
if(HAVE_ISO_FORTRAN_ENV)
  #add_compile_options(-DPSB_HAVE_ISO_FORTRAN_ENV)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DPSB_HAVE_ISO_FORTRAN_ENV")
  message(STATUS "-DPSB_HAVE_ISO_FORTRAN_ENV")
endif()
CHECK_Fortran_SOURCE_COMPILES(
  "flush(5); end"
  HAVE_FLUSH_STMT
  SRC_EXT f90
  )
if(HAVE_FLUSH_STMT)
  #add_compile_options(-DPSB_HAVE_FLUSH_STMT)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DPSB_HAVE_FLUSH_STMT")
  message(STATUS "-DPSB_HAVE_FLUSH_STMT")
endif()
CHECK_Fortran_SOURCE_COMPILES(
  "
module conftest_mod
  type foo
    integer :: i
  contains
    final  :: destroy_foo
  end type foo
  private destroy_foo
contains
  subroutine destroy_foo(a)
    type(foo) :: a
     ! Just a test
  end subroutine destroy_foo
end module conftest_mod
program conftest
  use conftest_mod
  type(foo) :: foovar
end program"
  HAVE_FINAL
  SRC_EXT f90
  )
if(HAVE_FINAL)
  # add_compile_options(-DPSB_HAVE_FINAL)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DPSB_HAVE_FINAL")
  message(STATUS "-DPSB_HAVE_FINAL")
endif()
CHECK_Fortran_SOURCE_COMPILES(
  "
program xtt
  type foo
    integer :: i
  end type foo
  type, extends(foo) :: new_foo
    integer :: j
  end type new_foo
  class(foo), allocatable  :: fooab
  type(new_foo) :: nfv
  integer :: info
  allocate(fooab, mold=nfv, stat=info)
end program"
  HAVE_MOLD
  SRC_EXT f90)
if(HAVE_MOLD)
  # add_compile_options(-DPSB_HAVE_MOLD)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DPSB_HAVE_MOLD")
  message(STATUS "-DPSB_HAVE_MOLD")
endif()
CHECK_Fortran_SOURCE_COMPILES(
  "
program conftest
  type foo
    integer :: i
  end type foo
  type, extends(foo) :: bar
    integer j
  end type bar
  type(bar) :: barvar
end program "
  HAVE_EXTENDS_TYPE_OF
  SRC_EXT f90)
if(HAVE_EXTENDS_TYPE_OF)
  # add_compile_options(-DPSB_HAVE_EXTENDS_TYPE_OF)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DPSB_HAVE_EXTENDS_TYPE_OF")
  message(STATUS "-DPSB_HAVE_EXTENDS_TYPE_OF")
endif()
CHECK_Fortran_SOURCE_COMPILES(
  "
program stt
  type foo
    integer :: i
  end type foo
  type, extends(foo) :: new_foo
    integer :: j
  end type new_foo
  type(foo) :: foov
  type(new_foo) :: nfv1, nfv2


  write(*,*) 'foov == nfv1? ', same_type_as(foov,nfv1)
  write(*,*) 'nfv2 == nfv1? ', same_type_as(nfv2,nfv1)
end program"
  HAVE_SAME_TYPE_AS
  SRC_EXT f90)
if(HAVE_SAME_TYPE_AS)
  # add_compile_options(-DPSB_HAVE_SAME_TYPE_AS)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DPSB_HAVE_SAME_TYPE_AS")
  message(STATUS "-DPSB_HAVE_SAME_TYPE_AS")
endif()

# ============================================================
# Robust MPI configuration (C + Fortran) for PSBLAS
# Works with modern CMake, MPICH/OpenMPI, local servers, and CI
# ============================================================

find_package(MPI REQUIRED COMPONENTS C Fortran)

if(MPI_FOUND)
  message(STATUS ">>> MPI found successfully!")
  message(STATUS "MPI version: ${MPI_VERSION}")

  # ------------------------------------------------------------
  # Extract include paths and library info from imported targets
  # ------------------------------------------------------------
  if(TARGET MPI::MPI_Fortran)
    get_target_property(_mpi_fortran_inc MPI::MPI_Fortran INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(_mpi_fortran_lib MPI::MPI_Fortran IMPORTED_LOCATION)
    get_target_property(_mpi_fortran_link MPI::MPI_Fortran INTERFACE_LINK_LIBRARIES)

    if(_mpi_fortran_inc)
      include_directories(BEFORE ${_mpi_fortran_inc})
      message(STATUS "MPI Fortran include paths: ${_mpi_fortran_inc}")
    endif()
  endif()

  if(TARGET MPI::MPI_C)
    get_target_property(_mpi_c_inc MPI::MPI_C INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(_mpi_c_lib MPI::MPI_C IMPORTED_LOCATION)

    if(_mpi_c_inc)
      include_directories(BEFORE ${_mpi_c_inc})
      message(STATUS "MPI C include paths: ${_mpi_c_inc}")
    endif()
  endif()

  # ------------------------------------------------------------
  # Fedora-specific workaround for noexecstack flag
  # ------------------------------------------------------------
  if((MPI_C_LINK_FLAGS MATCHES "noexecstack") OR (MPI_Fortran_LINK_FLAGS MATCHES "noexecstack"))
    message(WARNING
      "The `noexecstack` linker flag was found in MPI_<lang>_LINK_FLAGS.\n"
      "This can cause segmentation faults in Fortran codes.\n"
      "Replacing `noexecstack` with `execstack`."
    )
    string(REPLACE "noexecstack" "execstack" MPI_C_LINK_FLAGS "${MPI_C_LINK_FLAGS}")
    string(REPLACE "noexecstack" "execstack" MPI_Fortran_LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
  endif()

  # ------------------------------------------------------------
  # Compiler and linker flags
  # ------------------------------------------------------------
  set(CMAKE_C_COMPILE_FLAGS "${CMAKE_C_COMPILE_FLAGS} ${MPI_C_COMPILE_FLAGS}")
  set(CMAKE_C_LINK_FLAGS    "${CMAKE_C_LINK_FLAGS}    ${MPI_C_LINK_FLAGS}")
  set(CMAKE_Fortran_COMPILE_FLAGS "${CMAKE_Fortran_COMPILE_FLAGS} ${MPI_Fortran_COMPILE_FLAGS}")
  set(CMAKE_Fortran_LINK_FLAGS    "${CMAKE_Fortran_LINK_FLAGS}    ${MPI_Fortran_LINK_FLAGS}")

  message(STATUS "Fortran link flags: ${CMAKE_Fortran_LINK_FLAGS}")

  # ------------------------------------------------------------
  # Ensure mpi.mod is available for Fortran (legacy fallback)
  # ------------------------------------------------------------
  set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
  file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

  set(_mpi_mod_found FALSE)
  if(_mpi_fortran_inc)
    foreach(_mod_name mpi.mod MPI.mod)
      foreach(_inc ${_mpi_fortran_inc})
        if(EXISTS "${_inc}/${_mod_name}")
          file(COPY "${_inc}/${_mod_name}" DESTINATION "${CMAKE_Fortran_MODULE_DIRECTORY}")
          message(STATUS "Copied ${_mod_name} from ${_inc}")
          set(_mpi_mod_found TRUE)
          break()
        endif()
      endforeach()
      if(_mpi_mod_found)
        break()
      endif()
    endforeach()
  endif()

  if(NOT _mpi_mod_found)
    message(FATAL "mpi.mod not found in MPI include paths; assuming mpifort provides it internally. Skipping copy.")
  endif()


  # ------------------------------------------------------------
  # Enable MPI Fortran module support
  # ------------------------------------------------------------
  if(MPI_Fortran_HAVE_F90_MODULE OR MPI_Fortran_HAVE_F08_MODULE)
    add_compile_options(-DPSB_MPI_MOD)
    message(STATUS "Defined: -DPSB_MPI_MOD")
  endif()

  set(PSB_SERIAL_MPI OFF)

else()
  # ------------------------------------------------------------
  # Fallback to serial mode
  # ------------------------------------------------------------
  message(WARNING ">>> MPI not found — building in serial mode")
  add_compile_options(-DPSB_SERIAL_MPI -DPSB_MPI_MOD)
  set(PSB_SERIAL_MPI ON)
  set(CSERIALMPI "#define PSB_SERIAL_MPI")
endif()



#-------------------------------------------------------
# Find and Use OpenCoarrays IFF gfortran AND options set
#-------------------------------------------------------

if("${PSBLAS_USE_OpenCoarrays}" AND CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  message(STATUS "Set openCoarrays")
  find_package(OpenCoarrays)
endif()

#------------------------------
# Find Linear Algebra Libraries
#------------------------------
if(NOT APPLE)
  set(BLA_STATIC ON)
endif()
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
#add_compile_options(-DPSB_HAVE_LAPACK)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DPSB_HAVE_LAPACK")



#--------------------------------
# Find METIS partitioning library
#--------------------------------
include(${CMAKE_CURRENT_LIST_DIR}/cmake/FindMETIS.cmake)
find_package(METIS)


if(METIS_FOUND)
  message(STATUS "METIS PATH ${METIS_INCLUDES} and metis libraries ${METIS_LIBRARIES}")
  # Make sure this path is correct
#  set(METISINCFILE "metis.h")  # Adjust this to your actual path

  # Specify the configuration file
  #  set(HEADER_TEMPLATE "${CMAKE_CURRENT_SOURCE_DIR}/util/psb_metis_int.h.in")
  # set(HEADER_OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/include/psb_metis_int.h")

  # Configure the header file
  #configure_file(${HEADER_TEMPLATE} ${HEADER_OUTPUT} @ONLY)
 



  # Check for real sizes using try_compile
  include(CheckCSourceCompiles)

  # Function to check the size of a type
  function(check_metis_real_type type_name)
    set(source_code "
    #include <metis.h>
    #include <stdio.h>
    int main() {
        printf(\"%zu\\n\", sizeof(${type_name}));
        return 0;
    }")

    # Create a temporary source file
    file(WRITE "${CMAKE_BINARY_DIR}/CMakeFiles/CMakeTmp/test_size.c" "${source_code}")

    # Try to compile it
    try_compile(COMPILER_RESULT "${CMAKE_BINARY_DIR}/CMakeFiles/CMakeTmp"
        "${CMAKE_BINARY_DIR}/CMakeFiles/CMakeTmp/test_size.c")

    # Check the result and read the output
    if (COMPILER_RESULT)
        execute_process(COMMAND "${CMAKE_BINARY_DIR}/CMakeFiles/CMakeTmp/test_size"
            OUTPUT_VARIABLE type_size)
        string(STRIP "${type_size}" type_size)
        if (type_name STREQUAL "float")
            set(PSB_METIS_REAL_32 "${type_size}" PARENT_SCOPE)
	    #         add_definitions(-DPSB_METIS_REAL_32)
            set(CREALMETIS "#define PSB_METIS_REAL_32" PARENT_SCOPE)	
        elseif (type_name STREQUAL "double")
            set(PSB_METIS_REAL_64 "${type_size}" PARENT_SCOPE)
            #add_definitions(-DPSB_METIS_REAL_64)
            set(CREALMETIS "#define PSB_METIS_REAL_64" PARENT_SCOPE)		    
        endif()
    else()
        message(WARNING "Failed to compile test for type size: ${type_name}")
    endif()
  endfunction()

  # Check for both float and double
  check_metis_real_type(float)
  check_metis_real_type(double)

  # Set HAVE_METIS if METIS is found
  #add_compile_options(-DPSB_HAVE_METIS)
  #  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DPSB_HAVE_METIS")

  
  
  
  
  
  # Determine METIS_INDEX based on real type sizes
  if(DEFINED PSB_METIS_REAL_32)
    set(METIS_INDEX 32)
  elseif(DEFINED PSB_METIS_REAL_64)
    set(METIS_INDEX 64)
  else()
    message(WARNING "Neither METIS_REAL_32 nor METIS_REAL_64 is defined.")
    set(METIS_INDEX 64)  # Default to 64 if not defined
  endif()

   # Check conditions for LPK_SIZE and METIS_INDEX
  if(LPK_SIZE STREQUAL "4")
    if(METIS_INDEX STREQUAL "64")
        # Mismatch between METIS size and PSBLAS LPK
      message(FATAL " Mismatch between metis ${METIS_INDEX} size and psblas LPK size ${LPK_SIZE}")
    set(METIS_FOUND  FALSE)  
    endif()
  endif()

  if(LPK_SIZE STREQUAL "8")
    if(METIS_INDEX STREQUAL "32")
        # Mismatch between METIS size and PSBLAS LPK
            message(FATAL " Mismatch between metis ${METIS_INDEX} size and psblas LPK size ${LPK_SIZE}")
    set(METIS_FOUND  FALSE)
    endif()
  endif()

 if(METIS_FOUND)

    # Make sure this path is correct
  set(METISINCFILE "metis.h")  # Adjust this to your actual path

  # Specify the configuration file
  set(HEADER_TEMPLATE "${CMAKE_CURRENT_SOURCE_DIR}/util/psb_metis_int.h.in")
  set(HEADER_OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/include/psb_metis_int.h")

  # Configure the header file
  configure_file(${HEADER_TEMPLATE} ${HEADER_OUTPUT} @ONLY)

  # Set HAVE_METIS if METIS is found and coherent with the system settings
  #add_compile_options(-DPSB_HAVE_METIS)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DPSB_HAVE_METIS")


   
  set(CHAVEMETIS "#define PSB_HAVE_METIS")
  set(CINTMETIS "#define PSB_METIS_${METIS_INDEX}")
#  set(CREALMETIS "#define PSB_METIS_REAL_${LPK_SIZE}")
  # Configure the header file
  configure_file(${HEADER_TEMPLATE} ${HEADER_OUTPUT} @ONLY)
  # Set HAVE_METIS if METIS is found
  #add_compile_options(-DPSB_HAVE_METIS)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DPSB_HAVE_METIS")

endif()

endif()




#---------------------------------------------------
# Use standardized GNU install directory conventions
#---------------------------------------------------
include(GNUInstallDirs)
#set(mod_dir_tail "${${CMAKE_PROJECT_NAME}_dist_string}_${CMAKE_Fortran_COMPILER_ID}-${CMAKE_Fortran_COMPILER_VERSION}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/${${CMAKE_PROJECT_NAME}_dist_string}-tests")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
#set(CMAKE_INSTALL_LIBDIR "lib" CACHE STRING "Library install directory")
#set(CMAKE_INSTALL_INCLUDEDIR "include" CACHE STRING "Include directory")
#set(CMAKE_INSTALL_MODULDIR "modules" CACHE STRING "Modules directory")




#Ser variables exportable for other projects




message(STATUS "Initial CMAKE_INSTALL_LIBDIR: ${CMAKE_INSTALL_LIBDIR}")
set(PSB_CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})

if(NOT PSB_CMAKE_INSTALL_LIBDIR)
    message(STATUS "CMAKE_INSTALL_LIBDIR is set to default value lib")
    set(CMAKE_INSTALL_LIBDIR "lib" CACHE STRING "Library install directory" FORCE)
    set(PSB_CMAKE_INSTALL_LIBDIR ${CMAKE_INSTALL_LIBDIR})
else()
    set(CMAKE_INSTALL_LIBDIR ${PSB_CMAKE_INSTALL_LIBDIR})
    message(STATUS "CMAKE_INSTALL_LIBDIR is set to: ${CMAKE_INSTALL_LIBDIR}")
endif()

if(NOT PSB_CMAKE_INSTALL_INCLUDEDIR)
    message(STATUS "CMAKE_INSTALL_INCLUDEDIR is set to default value lib")
    set(CMAKE_INSTALL_INCLUDEDIR "include" CACHE STRING "Include directory" FORCE)
    set(PSB_CMAKE_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR})
else()
    set(CMAKE_INSTALL_INCLUDEDIR ${PSB_CMAKE_INSTALL_INCLUDEDIR})
    message(STATUS "CMAKE_INSTALL_INCLUDEDIR is set to: ${CMAKE_INSTALL_INCLUDEDIR}")
endif()

if(NOT PSB_CMAKE_INSTALL_MODULDIR)
    message(STATUS "CMAKE_INSTALL_MODULDIR is set to default value lib")
    set(CMAKE_INSTALL_MODULDIR "modules" CACHE STRING "Modules directory" FORCE)
    set(PSB_CMAKE_INSTALL_MODULDIR ${CMAKE_INSTALL_MODULDIR})
else()
    set(CMAKE_INSTALL_MODULDIR ${PSB_CMAKE_INSTALL_MODULDIR})
    message(STATUS "CMAKE_INSTALL_MODULDIR is set to: ${CMAKE_INSTALL_MODULDIR}")
endif()




#-----------------------------------
# Turn on testing/ctest capabilities
#-----------------------------------
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
#install(EXPORT ${CMAKE_PROJECT_NAME}-targets
#  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake"
#)

install(EXPORT ${CMAKE_PROJECT_NAME}-targets
  FILE ${CMAKE_PROJECT_NAME}Config.cmake
  NAMESPACE ${CMAKE_PROJECT_NAME}::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/"
)


include(CMakePackageConfigHelpers) # standard CMake module
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_PROJECT_NAME}ConfigVersion.cmake"
  VERSION "${psblas_VERSION}"
  COMPATIBILITY SameMajorVersion
  )

configure_file("${CMAKE_SOURCE_DIR}/cmake/${CMAKE_PROJECT_NAME}Config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${CMAKE_PROJECT_NAME}Config.cmake" @ONLY)

install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${CMAKE_PROJECT_NAME}Config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_PROJECT_NAME}ConfigVersion.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_PROJECT_NAME}Targets.cmake"
  DESTINATION
    "${CMAKE_INSTALL_LIBDIR}/cmake/${CMAKE_PROJECT_NAME}"
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

#----------------------------------
# Determine if we're using Open MPI
#---------------------------------
if(MPI_FOUND)
  execute_process(COMMAND ${MPIEXEC} --version
    OUTPUT_VARIABLE mpi_version_out)
  if (mpi_version_out MATCHES "[Oo]pen[ -][Mm][Pp][Ii]")
    message( STATUS "OpenMPI detected")
    set ( openmpi true )
  endif()
endif()

# Optionally check for CUDA requirement
option(PSB_BUILD_CUDA "Build CUDA code" OFF)


if(IPK_SIZE EQUAL 8)
  set(PSB_BUILD_CUDA OFF)
  message(STATUS "IPK8 is not compatible with CUDA. Cuda is now OFF ${PSB_BUILD_CUDA}") 
endif()


if(PSB_BUILD_CUDA)

    #if(NOT DEFINED PSB_CUDA_PATH)
    #        set(PSB_CUDA_PATH "/opt/cuda/12.8")
    #endif()
    # Include the CMakeLists for the cuda library
    include(${CMAKE_CURRENT_LIST_DIR}/cuda/CMakeLists.txt)
    include_directories("${PSB_CUDA_PATH}/include")
    message(STATUS "${PSB_CUDA_PATH}")


    #    find_package(CUDA REQUIRED)

    enable_language(CUDA)
    message(STATUS "Enabled CUDA ${CMAKE_CUDA_COMPILER_VERSION} ${CMAKE_CUDA_ARCHITECTURES};; ${CMAKE_CUDA_HOST_COMPILER_VERSION};")
    find_package(CUDAToolkit)
    message(STATUS "Enabled CUDA throguh find ${CUDAToolkit_VERSION_MAJOR} ${CUDAToolkit_VERSION};; ${CUDAToolkit_VERSION_MINOR};")
    #compute cuda versio for psblas
    math(EXPR PSB_CUDA_VERSION "${CUDAToolkit_VERSION_MAJOR} * 1000 + ${CUDAToolkit_VERSION_MINOR} * 10")



    

    message(STATUS "cuda version called has given ${PSB_CUDA_VERSION}:")

    # Check for CUDA version
    # set(PSB_CUDA_VERSION 12800)
    if(PSB_CUDA_VERSION)
            message(STATUS "CUDA version: ${PSB_CUDA_VERSION}")

        # Define macros for CUDA version
        #        add_definitions(-DPSB_HAVE_CUDA)
        #       add_definitions(-DPSB_CUDA_VERSION=${PSB_CUDA_VERSION})
        # math(EXPR PSB_CUDA_SHORT_VERSION "${PSB_CUDA_VERSION} / 1000")
        #        add_definitions(-DPSB_CUDA_SHORT_VERSION=${PSB_CUDA_SHORT_VERSION})
        set(PSB_CUDA_SHORT_VERSION ${CUDAToolkit_VERSION_MAJOR})
        message(STATUS "cuda version called has given ${PSB_CUDA_SHORT_VERSION}:")

        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DPSB_HAVE_CUDA")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DPSB_CUDA_VERSION=${PSB_CUDA_VERSION}")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DPSB_CUDA_SHORT_VERSION=${PSB_CUDA_SHORT_VERSION}")


        set(CHAVECUDA "#define PSB_HAVE_CUDA")
        set(CSHORTVCUDA "#define PSB_CUDA_SHORT_VERSION ${PSB_CUDA_SHORT_VERSION}")
        set(CVERSIONCUDA "#define PSB_CUDA_VERSION ${PSB_CUDA_VERSION}")

    else()
        message(FATAL_ERROR "CUDA version not found!")
    endif()
endif()

#------------------------------------------
# Configure the psb_config.h file
#------------------------------------------

message(STATUS "bin dir ${CMAKE_CURRENT_BINARY_DIR}; source dir ${CMAKE_CURRENT_SOURCE_DIR};;")
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/base/modules/psb_config.h.in 
    ${CMAKE_CURRENT_BINARY_DIR}/include/psb_config.h 	
    @ONLY  # Replace variables only
)



#---------------------------------------
# Add the PSBLAS libraries and utilities
#---------------------------------------

# Link order, left to right:
# cbind.a, util.a linsolve.a prec.a base.a





include(${CMAKE_CURRENT_LIST_DIR}/base/CMakeLists.txt)

include_directories("${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}")


  foreach(path IN LISTS base_header_C_files)

        # Copy the header file to the include directory
        file(COPY "${path}" DESTINATION "${CMAKE_BINARY_DIR}/include")
        

  endforeach()

if(WIN32)
  add_library(psb_base_C STATIC ${base_source_C_files})
  target_compile_definitions(psb_base_C
    PRIVATE -DWIN32 -D_LIB -DWIN64)
  set_target_properties(psb_base_C
    PROPERTIES
    LINKER_LANGUAGE C
    POSITION_INDEPENDENT_CODE TRUE)
  target_link_libraries(psb_base_C
    PUBLIC kernel32 user32 shell32)
  add_library(base ${base_source_files})
  target_link_libraries(base
    PUBLIC psb_base_C)
else()
  add_library(base_C OBJECT ${base_source_C_files})
  add_library(base ${base_source_files} $<TARGET_OBJECTS:base_C>)
endif()


# Set the Fortran module output directory for all targets
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
#set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/include")



message(STATUS "fortran module direcotry ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}")


include_directories(${MPI_Fortran_INCLUDE_PATH})


message(STATUS "Using MPI include at: ${MPI_Fortran_INCLUDE_PATH}")
    


set_target_properties(base
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules"  
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME psb_base
  LINKER_LANGUAGE Fortran
  )

target_include_directories(base PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/modules> 
  $<INSTALL_INTERFACE:modules>)

message(STATUS "include dir := ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}")

#target_include_directories(base PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})


target_link_libraries(base
  PUBLIC ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES} ${LAPACK95_LIBRARIES}
  PUBLIC ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${BLAS95_LIBRARIES})

#add_custom_command(
#    TARGET base POST_BUILD
#    COMMAND ${CMAKE_COMMAND} -E cmake_copy_f90_mod
#    ${MPI_Fortran_INCLUDE_PATH}mpi.mod ${CMAKE_BINARY_DIR}/include/mpi.mod
#)




include(${CMAKE_CURRENT_LIST_DIR}/prec/CMakeLists.txt)
add_library(prec ${prec_source_files})
set_target_properties(prec
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules"
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME psb_prec
  LINKER_LANGUAGE Fortran
  )
target_include_directories(prec PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/modules>
  $<INSTALL_INTERFACE:modules>)
target_link_libraries(prec PUBLIC base)


include(${CMAKE_CURRENT_LIST_DIR}/linsolve/CMakeLists.txt)
add_library(linsolve ${linsolve_source_files})
set_target_properties(linsolve
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules"
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME psb_linsolve
  LINKER_LANGUAGE Fortran
  )
target_include_directories(linsolve PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/modules>
  $<INSTALL_INTERFACE:modules>)
target_link_libraries(linsolve PUBLIC base prec)



include(${CMAKE_CURRENT_LIST_DIR}/ext/CMakeLists.txt)
add_library(ext ${ext_source_files})
set_target_properties(ext
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules"
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME psb_ext
  LINKER_LANGUAGE Fortran
  )
target_include_directories(ext PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/modules>
  $<INSTALL_INTERFACE:modules>)
target_link_libraries(ext PUBLIC base prec)  #TODO: check actual dependencies




include(${CMAKE_CURRENT_LIST_DIR}/util/CMakeLists.txt)

if(WIN32)
  if(METIS_FOUND)
    add_library(psb_util_C STATIC ${util_source_C_files})
    target_compile_definitions(psb_util_C
      PRIVATE -DWIN32 -D_LIB -DWIN64)
    set_target_properties(psb_util_C
      PROPERTIES
      LINKER_LANGUAGE C
      POSITION_INDEPENDENT_CODE TRUE)
    target_link_libraries(psb_util_C
      PUBLIC kernel32 user32 shell32)
  endif()
  add_library(util ${util_source_files})
  if(METIS_FOUND)
    target_link_libraries(util
      PUBLIC psb_util_C)
  endif()
else()

  if(METIS_FOUND)
    foreach(file IN LISTS util_source_C_metis_files)
      list(APPEND util_source_C_files file)
    endforeach()
  endif()
  add_library(psb_util_C OBJECT ${util_source_C_files})

  add_library(util ${util_source_files} $<TARGET_OBJECTS:psb_util_C>)
endif()
set_target_properties(util
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules"
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME psb_util
  LINKER_LANGUAGE Fortran
  )
target_include_directories(util PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/modules>
  $<INSTALL_INTERFACE:modules>)
target_link_libraries(util PUBLIC base prec)

if(METIS_FOUND)

  target_include_directories(util
    PUBLIC ${METIS_INCLUDES})
  target_include_directories(psb_util_C
    PUBLIC ${METIS_INCLUDES})
  target_link_libraries(util
    PUBLIC ${METIS_LIBRARIES})
#  target_compile_definitions(psb_util_C
#    PUBLIC PSB_HAVE_METIS) #TDDO: CHECK IF THAT _ IS CORRECT
#  target_compile_definitions(util
#    PUBLIC PSB_HAVE_METIS)
endif()









# Include headers from the 'include' directory in the current directory
include_directories(${CMAKE_BINARY_DIR}/include)


include(${CMAKE_CURRENT_LIST_DIR}/cbind/CMakeLists.txt)
if(WIN32)
  add_library(psb_cbind_C STATIC ${cbind_source_C_files})
  target_compile_definitions(psb_cbind_C
    PRIVATE -DWIN32 -D_LIB -DWIN64)
  set_target_properties(psb_cbind_C
    PROPERTIES
    LINKER_LANGUAGE C
    POSITION_INDEPENDENT_CODE TRUE)
  target_link_libraries(psb_cbind_C
    PUBLIC kernel32 user32 shell32)
  add_library(cbind ${cbind_source_files})
  target_link_libraries(cbind
    PUBLIC psb_cbind_C)
else()
  add_library(cbind_C OBJECT ${cbind_source_C_files})
  add_library(cbind ${cbind_source_files} $<TARGET_OBJECTS:cbind_C>)
  
endif()


#add_library(cbind ${cbind_source_files})
set_target_properties(cbind
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules"
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME psb_cbind
  LINKER_LANGUAGE Fortran
  )
#target_include_directories(cbind PUBLIC
#  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/modules>
#  $<INSTALL_INTERFACE:modules>)
# Include directories for the cbind library
target_include_directories(cbind PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/modules>  # Path for building
  $<INSTALL_INTERFACE:modules>  # Path for installation
)

target_link_libraries(cbind PUBLIC base prec linsolve ext util)



# Custom command to copy all header files
#add_custom_command(
#  OUTPUT ${CMAKE_BINARY_DIR}/include/  # Dummy output to represent the target directory
#  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/include  # Create the include directory if it doesn't exist
#  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${cbind_header_C_files} ${CMAKE_BINARY_DIR}/include/  # Copy all headers
#  DEPENDS ${cbind_header_C_files}  # Make the copy depend on the header files
#  COMMENT "Copying header files to include directory"
#)

# Create a custom target to copy headers
#add_custom_target(copy_headers ALL DEPENDS ${CMAKE_BINARY_DIR}/include/)

  foreach(path IN LISTS cbind_header_C_files)

        # Copy the header file to the include directory
        file(COPY "${path}" DESTINATION "${CMAKE_BINARY_DIR}/include")
        

  endforeach()
message(STATUS "Copied .h files to ${CMAKE_BINARY_DIR}/include")


#target_include_directories(cbind PUBLIC
#  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/modules>
#  $<INSTALL_INTERFACE:modules>)
# Include directories for the cbind library
#target_include_directories(cbind_C PUBLICF
#  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>  # Path for building
#  $<INSTALL_INTERFACE:include>  # Path for installation
#  ${CMAKE_BINARY_DIR}/include  # Include the copied headers
#)






#########################################
#######	BUILD CUDA LIBRARY ##############
#########################################


# Optionally check for CUDA requirement
#option(PSB_BUILD_CUDA "Build CUDA code" OFF)

if(PSB_BUILD_CUDA)

	#    if(NOT DEFINED PSB_CUDA_PATH)
	#   set(PSB_CUDA_PATH "/opt/cuda/12.8")    
	#endif()
    # Include the CMakeLists for the cbind library
    #include(${CMAKE_CURRENT_LIST_DIR}/cuda/CMakeLists.txt)
    #include_directories("${PSB_CUDA_PATH}/include")
    #message(STATUS "${PSB_CUDA_PATH}")


    #    find_package(CUDA REQUIRED)

    #enable_language(CUDA)
    
    # Check for CUDA version
    #set(PSB_CUDA_VERSION 12800)
    #if(PSB_CUDA_VERSION)
    #  message(STATUS "CUDA version: ${PSB_CUDA_VERSION}")

	# Define macros for CUDA version
	#        add_definitions(-DPSB_HAVE_CUDA)
	#	add_definitions(-DPSB_CUDA_VERSION=${PSB_CUDA_VERSION})
	#math(EXPR PSB_CUDA_SHORT_VERSION "${PSB_CUDA_VERSION} / 1000")
	#        add_definitions(-DPSB_CUDA_SHORT_VERSION=${PSB_CUDA_SHORT_VERSION})
       
	#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DPSB_HAVE_CUDA")
	#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DPSB_CUDA_VERSION=${PSB_CUDA_VERSION}")
	#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DPSB_CUDA_SHORT_VERSION=${PSB_CUDA_SHORT_VERSION}")


	#set(CHAVECUDA "#define PSB_HAVE_CUDA")
	#set(CSHORTVCUDA "#define PSB_CUDA_SHORT_VERSION ${PSB_CUDA_SHORT_VERSION}")
	#set(CVERSIONCUDA "#define PSB_CUDA_VERSION ${PSB_CUDA_VERSION}")

	#else()
	#message(FATAL_ERROR "CUDA version not found!")
	#endif()



    # Define the CUDA library
    #if(WIN32)
        #add_library(psb_cuda_C STATIC ${cuda_source_files})
        #target_compile_definitions(psb_cuda_C
        #    PRIVATE -DWIN32 -D_LIB -DWIN64)
        #set_target_properties(psb_cuda_C
        #    PROPERTIES
        #    LINKER_LANGUAGE C
        #    POSITION_INDEPENDENT_CODE TRUE)
        #target_link_libraries(psb_cuda_C
        #    PUBLIC kernel32 user32 shell32)
    #else()
        #add_library(psb_cuda_C OBJECT ${cuda_source_files})
    #endif()
   
    foreach(path IN LISTS cuda_header_C_files)

        # Copy the header file to the include directory
        file(COPY "${path}" DESTINATION "${CMAKE_BINARY_DIR}/include")
        

    endforeach()
    message(STATUS "Copied .h files to ${CMAKE_BINARY_DIR}/include")
    foreach(path IN LISTS cuda_header_cu_files)

        # Copy the header file to the include directory
        file(COPY "${path}" DESTINATION "${CMAKE_BINARY_DIR}/include")


    endforeach()
    message(STATUS "Copied .cuh files to ${CMAKE_BINARY_DIR}/include")




    add_library(psb_cuda_C OBJECT ${cuda_source_C_files} ${cuda_source_cu_files})
    

    # Create the main CUDA library
    add_library(cuda ${cuda_source_files})

    # Set properties for the CUDA library
    set_target_properties(cuda
        PROPERTIES
        POSITION_INDEPENDENT_CODE TRUE
        OUTPUT_NAME psb_cuda
	LINKER_LANGUAGE C)
    
    # Include directories for the CUDA library
    target_include_directories(cuda PUBLIC
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/modules>  # Path for building
        $<INSTALL_INTERFACE:modules>  # Path for installation
        #/opt/cuda/12.8/include
    )

    # Link with other necessary libraries
    target_link_libraries(cuda PUBLIC base prec linsolve ext util)
endif()




#if(MPI_FOUND)
# Copy mpi.mod from the first available path in MPI_Fortran_INCLUDE_PATH
#  set(MPI_MOD_COPIED FALSE)

#  foreach(path IN LISTS MPI_Fortran_INCLUDE_PATH)
    # Construct the full path to the mpi.mod file
#    set(mpi_mod_path "${path}/mpi.mod")
    
    # Check if the mpi.mod file exists
#    if(EXISTS "${mpi_mod_path}")
        # Copy the mpi.mod file to the module directory
#        file(COPY "${mpi_mod_path}" DESTINATION "${CMAKE_Fortran_MODULE_DIRECTORY}")
        #message(STATUS "Copied mpi.mod from ${mpi_mod_path} to ${CMAKE_Fortran_MODULE_DIRECTORY}")
        #set(MPI_MOD_COPIED TRUE)
        #break()  # Exit the loop once we've copied the file
    #endif()
  #endforeach()

  #if(NOT MPI_MOD_COPIED)
  #  message(WARNING "mpi.mod not found in any of the specified paths: ${MPI_Fortran_INCLUDE_PATH}")
  #endif()
    
    

  #foreach(lib base prec linsolve ext util cbind)

   # target_link_libraries(${lib} PUBLIC ${MPI_C_LIBRARIES} ${MPI_Fortran_LIBRARIES})
  #endforeach()
#endif()

if(OpenCoarrays_FOUND)
  foreach(lib base prec linsolve ext util cbind)  #TODO: check if cbind goes here!
    target_link_libraries(${lib} PUBLIC OpenCoarrays::caf_mpi_static)
  endforeach()
endif()





message(STATUS "\t\t ${CMAKE_INSTALL_LIBDIR}")
# Install the header files to the include directory
#install(FILES ${cbind_header_C_files}
#  DESTINATION include
#)

#install(DIRECTORY "${CMAKE_BINARY_DIR}/include" DESTINATION "include"
#  FILES_MATCHING PATTERN "*.h")

#install(DIRECTORY "${CMAKE_BINARY_DIR}/modules" DESTINATION "modules"
#  FILES_MATCHING PATTERN "*.mod")


#install(DIRECTORY "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}" DESTINATION "include"
#  FILES_MATCHING PATTERN "*.h")

#install(DIRECTORY "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_MODULDIR}" DESTINATION "modules"
#  FILES_MATCHING PATTERN "*.mod")

# Install header files
install(DIRECTORY ${CMAKE_BINARY_DIR}/include/
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"  # This will place headers in <prefix>/include
  FILES_MATCHING PATTERN "*.h"
)

# Install module files
install(DIRECTORY ${CMAKE_BINARY_DIR}/modules/
  DESTINATION "${CMAKE_INSTALL_MODULDIR}"  # This will place .mod files in <prefix>/modules
  FILES_MATCHING PATTERN "*.mod"
)



install(TARGETS base prec linsolve ext util cbind
  EXPORT ${CMAKE_PROJECT_NAME}-targets
  DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  )
 
if(PSB_BUILD_CUDA)
 install(TARGETS cuda
   EXPORT ${CMAKE_PROJECT_NAME}-targets
   DESTINATION "${CMAKE_INSTALL_LIBDIR}"
   LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
   )

endif() 

if(WIN32)
  install(TARGETS psb_base_C
    EXPORT ${CMAKE_PROJECT_NAME}-targets
    DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    )
  if(METIS_FOUND)
    install(TARGETS psb_util_C
      EXPORT ${CMAKE_PROJECT_NAME}-targets
      DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      )
  endif()
endif()


#configure_package_config_file(
#  INSTALL_DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/psblas"
#)

#install(FILES 
#  "${CMAKE_CURRENT_BINARY_DIR}/psblasConfig.cmake"
#  "${CMAKE_CURRENT_BINARY_DIR}/psblasConfigVersion.cmake"
#  DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake/psblas"
#)


export(
  EXPORT ${CMAKE_PROJECT_NAME}-targets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/psblasTargets.cmake"
    NAMESPACE ${CMAKE_PROJECT_NAME}::
)




# Set the installation directory for the test files
set(INSTALL_TEST_DIR "${CMAKE_INSTALL_PREFIX}/samples" CACHE PATH "Installation directory for sample files")

function(install_directory_recursive source_dir install_base_dir) # Function to install a directory and its subdirectories recursively
    file(GLOB_RECURSE ALL_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/${source_dir}" "${source_dir}/*")

    foreach(FILE_PATH IN LISTS ALL_FILES)
        # Construct the full source and destination paths
        set(FULL_SOURCE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/${source_dir}/${FILE_PATH}")
        set(FULL_INSTALL_PATH "${install_base_dir}/${FILE_PATH}")

        # Check if it's a directory
        if(IS_DIRECTORY "${FULL_SOURCE_PATH}")
            # Create the directory in the install destination
            file(MAKE_DIRECTORY "${FULL_INSTALL_PATH}")
        else()
            # Install the file
            install(FILES "${FULL_SOURCE_PATH}" DESTINATION "${install_base_dir}" RENAME "${FILE_PATH}")
        endif()
    endforeach()
endfunction()

# Install cbind/test directory
install_directory_recursive(cbind/test "${INSTALL_TEST_DIR}/cbind")

# Install test/fileread directory
install_directory_recursive(test/fileread "${INSTALL_TEST_DIR}/fileread")

# Install test/pdegen directory
install_directory_recursive(test/pdegen "${INSTALL_TEST_DIR}/pdegen")


message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX} - ${PSB_CMAKE_INSTALL_PREFIX};")
message(STATUS "CMAKE_INSTALL_LIBDIR: ${CMAKE_INSTALL_LIBDIR} - ${PSB_CMAKE_INSTALL_LIBDIR};")
message(STATUS "CMAKE_INSTALL_INCLUDEDIR: ${CMAKE_INSTALL_INCLUDEDIR} - ${PSB_CMAKE_INSTALL_INCLUDEDIR};")
message(STATUS "CMAKE_INSTALL_MODULDIR: ${CMAKE_INSTALL_MODULDIR} - ${PSB_CMAKE_INSTALL_MODULDIR};")




#-----------------
# Add PSBLAS tests
#-----------------

# Unit tests targeting each function, argument, and branch of code
# add_mpi_test(initialize_mpi 2 initialize_mpi)

