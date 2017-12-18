if (NOT MPI_C_FOUND)
  find_package(MPI REQUIRED)

  set(CMAKE_C_COMPILE_FLAGS ${CMAKE_C_COMPILE_FLAGS} ${MPI_C_COMPILE_FLAGS})
  set(CMAKE_C_LINK_FLAGS ${CMAKE_C_LINK_FLAGS} ${MPI_C_LINK_FLAGS})
  set(CMAKE_Fortran_COMPILE_FLAGS ${CMAKE_Fortran_COMPILE_FLAGS} ${MPI_Fortran_COMPILE_FLAGS})
  set(CMAKE_Fortran_LINK_FLAGS ${CMAKE_Fortran_LINK_FLAGS} ${MPI_Fortran_LINK_FLAGS})
  include_directories(BEFORE ${MPI_C_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH})
endif()

if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  set(gfortran_compiler true)
endif()

include(CheckIncludeFile)
CHECK_INCLUDE_FILE("alloca.h" HAVE_ALLOCA)
if(NOT HAVE_ALLOCA)
  add_definitions(-DALLOCA_MISSING)
  message(WARNING "Could not find <alloca.h>. Assuming functionality is provided elsewhere.")
endif()

set(old_cmake_required_includes "${CMAKE_REQUIRED_INCLUDES}")
if(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES};${MPI_C_INCLUDE_PATH})
else()
  set(CMAKE_REQUIRED_INCLUDES ${MPI_C_INCLUDE_PATH})
endif()
set(old_cmake_required_flags "${CMAKE_REQUIRED_FLAGS}")
if(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS};${MPI_C_COMPILE_FLAGS};${MPI_C_LINK_FLAGS})
else()
  set(CMAKE_REQUIRED_FLAGS ${MPI_C_COMPILE_FLAGS};${MPI_C_LINK_FLAGS})
endif()
set(old_cmake_required_libraries "${CMAKE_REQUIRED_LIBRARIES}")
if(CMAKE_REQUIRED_LIBRARIES)
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES};${MPI_C_LIBRARIES})
else()
  set(CMAKE_REQUIRED_LIBRARIES ${MPI_C_LIBRARIES})
endif()

set(MPI_HEADERS mpi.h)
include(CheckIncludeFiles)
CHECK_INCLUDE_FILES("mpi.h;mpi-ext.h" HAVE_MPI_EXT)
if(HAVE_MPI_EXT)
  add_definitions(-DHAVE_MPI_EXT_H)
  set(MPI_HEADERS ${MPI_HEADERS};mpi-ext.h)
endif()

set(CMAKE_REQUIRED_INCLUDES ${old_cmake_required_includes})
set(CMAKE_REQUIRED_FLAGS ${old_cmake_required_flags})
set(CMAKE_REQUIRED_LIBRARIES ${old_cmake_required_libraries})

#---------------------------------------------------
# Windows Intel MPI compatibility, see GH issue #435
#---------------------------------------------------
set(old_cmake_required_includes "${CMAKE_REQUIRED_INCLUDES}")
if(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES};${MPI_C_INCLUDE_PATH})
else()
  set(CMAKE_REQUIRED_INCLUDES ${MPI_C_INCLUDE_PATH})
endif()
CHECK_INCLUDE_FILES("mpi.h" HAVE_MPI_H)

include(CheckSymbolExists)
CHECK_SYMBOL_EXISTS(I_MPI_VERSION "mpi.h" HAVE_Intel_MPI)
if(HAVE_Intel_MPI AND WIN32)
  add_definitions(-DUSE_GCC)
endif()
set(CMAKE_REQUIRED_INCLUDES ${old_cmake_required_includes})

set(PSBLAS_SO_VERSION 0)
if(gfortran_compiler)
  if(NOT CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 7.0.0)
    set(PSBLAS_SO_VERSION 2)
  elseif(NOT CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 6.0.0)
    set(PSBLAS_SO_VERSION 1)
  endif()
endif()

set(PSBLAS_MODDIR "${CMAKE_INSTALL_INCLUDEDIR}/${mod_dir_tail}")
set(MOD_DIR_FLAG "${CMAKE_Fortran_MODDIR_FLAG}")
set(MPI_LIBS "")
foreach( lib IN LISTS MPI_Fortran_LIBRARIES)
  set(MPI_LIBS "${MPI_LIBS} \"${lib}\"")
endforeach()
string(STRIP "${MPI_LIBS}" MPI_LIBS)

set(prec_source_files
./psb_prec_const_mod.f90
./psb_c_prec_mod.f90
./psb_c_base_prec_mod.f90
./psb_d_bjacprec.f90
./psb_z_base_prec_mod.f90
./psb_d_diagprec.f90
./impl/psb_dprecbld.f90
./impl/psb_cprecinit.f90
./impl/psb_dprecinit.f90
./impl/psb_zilu_fct.f90
./impl/psb_silu_fct.f90
./impl/psb_cprecset.f90
./impl/psb_s_prec_type_impl.f90
./impl/psb_dilu_fct.f90
./impl/psb_cprecbld.f90
./impl/psb_c_diagprec_impl.f90
./impl/psb_s_bjacprec_impl.f90
./impl/psb_z_prec_type_impl.f90
./impl/psb_d_nullprec_impl.f90
./impl/psb_c_nullprec_impl.f90
./impl/psb_d_diagprec_impl.f90
./impl/psb_zprecinit.f90
./impl/psb_d_prec_type_impl.f90
./impl/psb_c_prec_type_impl.f90
./impl/psb_s_nullprec_impl.f90
./impl/psb_sprecinit.f90
./impl/psb_z_bjacprec_impl.f90
./impl/psb_d_bjacprec_impl.f90
./impl/psb_z_nullprec_impl.f90
./impl/psb_dprecset.f90
./impl/psb_c_bjacprec_impl.f90
./impl/psb_zprecset.f90
./impl/psb_sprecbld.f90
./impl/psb_cilu_fct.f90
./impl/psb_z_diagprec_impl.f90
./impl/psb_sprecset.f90
./impl/psb_zprecbld.f90
./impl/psb_s_diagprec_impl.f90
./psb_s_base_prec_mod.f90
./psb_s_bjacprec.f90
./psb_c_prec_type.f90
./psb_s_nullprec.f90
./psb_d_nullprec.f90
./psb_prec_mod.f90
./psb_d_prec_type.f90
./psb_z_prec_type.f90
./psb_z_prec_mod.f90
./psb_z_nullprec.f90
./psb_c_diagprec.f90
./psb_z_bjacprec.f90
./psb_s_prec_type.f90
./psb_c_nullprec.f90
./psb_prec_type.f90
./psb_d_prec_mod.f90
./psb_s_diagprec.f90
./psb_c_bjacprec.f90
./psb_s_prec_mod.f90
./psb_z_diagprec.f90
./psb_d_base_prec_mod.f90
)

add_library(prec ${prec_source_files})


target_link_libraries(prec PRIVATE ${MPI_C_LIBRARIES} ${MPI_Fortran_LIBRARIES} base)

install(DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  FILES_MATCHING PATTERN "*.mod")

set_target_properties ( prec
  PROPERTIES
  SOVERSION ${PSBLAS_SO_VERSION}
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${mod_dir_tail}"
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME psb_prec
)

install(TARGETS prec EXPORT PSBLAS_Targets
  DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)
