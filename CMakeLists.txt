# FDEFINES : -DHAVE_LAPACK -DHAVE_FINAL -DHAVE_ISO_FORTRAN_ENV -DHAVE_FLUSH_STMT -DHAVE_VOLATILE -DSERIAL_MPI -DMPI_MOD
# CDEFINES : -DLowerUnderscore -DPtr64Bits

#-----------------------------------
# Set oldest allowable CMake version
#-----------------------------------
cmake_minimum_required(VERSION 3.10)
cmake_policy(VERSION 3.10.0...3.12.4)

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
message( STATUS "Building PSBLAS version: ${full_git_describe}" )

#------------------------------------------
# Name project and specify source languages
#------------------------------------------
project(psblas
  VERSION "${PSBLAS_Version}"
  LANGUAGES C Fortran)

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

#------------------------------------
# Fortran name mangling introspection
#------------------------------------
include(FortranCInterface)
string(TOUPPER ${FortranCInterface_GLOBAL__CASE} fc_case)
message(STATUS "Name mangling capitalization: ${fc_case}")
message(STATUS "Name mangling fortran global suffix underscore: ${FortranCInterface_GLOBAL__SUFFIX}")
if(FortranCInterface_GLOBAL__SUFFIX STREQUAL "")
  add_definitions("-D${fc_case}Case")
elseif(FortranCInterface_GLOBAL__SUFFIX STREQUAL "_")
  add_definitions("-D${fc_case}Underscore")
elseif(FortranCInterface_GLOBAL__SUFFIX STREQUAL "__")
  add_definitions("-D${fc_case}DoubleUnderscore")
else()
  message( FATAL_ERROR "Fortran name mangling suffix, \'${FortranCInterface_GLOBAL__SUFFIX}\', unknown to PSBLAS")
endif()

if(NOT ${WIN32})
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
endif()

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
  add_definitions(-DHAVE_FINAL)
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
  add_definitions(-DHAVE_MOLD)
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
  add_definitions(-DHAVE_EXTENDS_TYPE_OF)
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
  add_definitions(-DHAVE_SAME_TYPE_AS)
endif()

#----------------------------------------------------------------------------
# Find MPI and set some flags so that FC and CC can point to gfortran and gcc
#----------------------------------------------------------------------------
find_package( MPI )

if(MPI_FOUND)
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

  #----------------
  # Setup MPI flags
  #----------------
  list(REMOVE_DUPLICATES MPI_Fortran_INCLUDE_PATH)
  set(CMAKE_C_COMPILE_FLAGS ${CMAKE_C_COMPILE_FLAGS} ${MPI_C_COMPILE_FLAGS})
  set(CMAKE_C_LINK_FLAGS ${CMAKE_C_LINK_FLAGS} ${MPI_C_LINK_FLAGS})
  set(CMAKE_Fortran_COMPILE_FLAGS ${CMAKE_Fortran_COMPILE_FLAGS} ${MPI_Fortran_COMPILE_FLAGS})
  set(CMAKE_Fortran_LINK_FLAGS ${CMAKE_Fortran_LINK_FLAGS} ${MPI_Fortran_LINK_FLAGS})
  include_directories(BEFORE ${MPI_C_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH})

  if(MPI_Fortran_HAVE_F90_MODULE OR MPI_Fortran_HAVE_F08_MODULE)
    add_definitions(-DMPI_MOD)
  endif()
  set(SERIAL_MPI OFF)
else()
  add_definitions(-DSERIAL_MPI)
  add_definitions(-DMPI_MOD)
  set(SERIAL_MPI ON)
endif()

#------------------------------
# Find Linear Algebra Libraries
#------------------------------
#set(BLA_STATIC 1)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
add_definitions(-DHAVE_LAPACK)


#--------------------------------
# Find METIS partitioning library
#--------------------------------
include(${CMAKE_CURRENT_LIST_DIR}/cmake/FindMETIS.cmake)
find_package(METIS)

#---------------------------------------------------
# Use standardized GNU install directory conventions
#---------------------------------------------------
include(GNUInstallDirs)
#set(mod_dir_tail "${${CMAKE_PROJECT_NAME}_dist_string}_${CMAKE_Fortran_COMPILER_ID}-${CMAKE_Fortran_COMPILER_VERSION}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/${${CMAKE_PROJECT_NAME}_dist_string}-tests")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

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


if(MPI_FOUND)
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
  execute_process(COMMAND ${MPIEXEC} --version
    OUTPUT_VARIABLE mpi_version_out)
  if (mpi_version_out MATCHES "[Oo]pen[ -][Mm][Pp][Ii]")
    message( STATUS "OpenMPI detected")
    set ( openmpi true )
  endif()
endif()

#---------------------------------------
# Add the PSBLAS libraries and utilities
#---------------------------------------

# Link order, left to right:
# cbind.a, util.a krylov.a prec.a base.a

include(${CMAKE_CURRENT_LIST_DIR}/base/CMakeLists.txt)
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
set_target_properties(base
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}"
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME psb_base
  LINKER_LANGUAGE Fortran
  )
target_include_directories(base PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
if(MPI_FOUND)
  target_link_libraries(base PUBLIC ${MPI_C_LIBRARIES} ${MPI_Fortran_LIBRARIES})
endif()
target_link_libraries(base
  PUBLIC ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES} ${LAPACK95_LIBRARIES}
  PUBLIC ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${BLAS95_LIBRARIES})

include(${CMAKE_CURRENT_LIST_DIR}/prec/CMakeLists.txt)
add_library(prec ${prec_source_files})
set_target_properties(prec
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}"
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME psb_prec
  LINKER_LANGUAGE Fortran
  )
target_include_directories(prec PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(prec PUBLIC base)

include(${CMAKE_CURRENT_LIST_DIR}/krylov/CMakeLists.txt)
add_library(krylov ${krylov_source_files})
set_target_properties(krylov
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}"
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME psb_krylov
  LINKER_LANGUAGE Fortran
  )
target_include_directories(krylov PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(krylov PUBLIC base prec)

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
  add_library(psb_util_C OBJECT ${util_source_C_files})
  add_library(util ${util_source_files} $<TARGET_OBJECTS:psb_util_C>)
endif()
set_target_properties(util
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}"
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME psb_util
  LINKER_LANGUAGE Fortran
  )
target_include_directories(util PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(util PUBLIC base prec)
if(METIS_FOUND)
  message(STATUS ${METIS_INCLUDES})
  target_include_directories(util
    PUBLIC ${METIS_INCLUDES})
  target_include_directories(psb_util_C
    PUBLIC ${METIS_INCLUDES})
  target_link_libraries(util
    PUBLIC ${METIS_LIBRARIES})
  target_compile_definitions(psb_util_C
    PUBLIC HAVE_METIS_)
  target_compile_definitions(util
    PUBLIC HAVE_METIS)
endif()

install(DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  FILES_MATCHING PATTERN "*.mod")
install(TARGETS base prec krylov util
  EXPORT ${CMAKE_PROJECT_NAME}-targets
  DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  )
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

#-----------------
# Add PSBLAS tests
#-----------------

# Unit tests targeting each function, argument, and branch of code
# add_mpi_test(initialize_mpi 2 initialize_mpi)
