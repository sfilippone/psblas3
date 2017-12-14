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

add_library(util SHARED
metis_int.c
psb_amd_order.c
psb_blockpart_mod.f90
psb_c_hbio_impl.f90
psb_c_mat_dist_impl.f90
psb_c_mat_dist_mod.f90
psb_c_mmio_impl.f90
psb_c_renum_impl.F90
psb_d_hbio_impl.f90
psb_d_mat_dist_impl.f90
psb_d_mat_dist_mod.f90
psb_d_mmio_impl.f90
psb_d_renum_impl.F90
psb_gps_mod.f90
psb_hbio_mod.f90
psb_i_mmio_impl.f90
psb_mat_dist_mod.f90
psb_metispart_mod.F90
psb_mmio_mod.F90
psb_renum_mod.f90
psb_s_hbio_impl.f90
psb_s_mat_dist_impl.f90
psb_s_mat_dist_mod.f90
psb_s_mmio_impl.f90
psb_s_renum_impl.F90
psb_util_mod.f90
psb_z_hbio_impl.f90
psb_z_mat_dist_impl.f90
psb_z_mat_dist_mod.f90
psb_z_mmio_impl.f90
psb_z_renum_impl.F90
)

add_library(util_static STATIC
metis_int.c
psb_amd_order.c
psb_blockpart_mod.f90
psb_c_hbio_impl.f90
psb_c_mat_dist_impl.f90
psb_c_mat_dist_mod.f90
psb_c_mmio_impl.f90
psb_c_renum_impl.F90
psb_d_hbio_impl.f90
psb_d_mat_dist_impl.f90
psb_d_mat_dist_mod.f90
psb_d_mmio_impl.f90
psb_d_renum_impl.F90
psb_gps_mod.f90
psb_hbio_mod.f90
psb_i_mmio_impl.f90
psb_mat_dist_mod.f90
psb_metispart_mod.F90
psb_mmio_mod.F90
psb_renum_mod.f90
psb_s_hbio_impl.f90
psb_s_mat_dist_impl.f90
psb_s_mat_dist_mod.f90
psb_s_mmio_impl.f90
psb_s_renum_impl.F90
psb_util_mod.f90
psb_z_hbio_impl.f90
psb_z_mat_dist_impl.f90
psb_z_mat_dist_mod.f90
psb_z_mmio_impl.f90
psb_z_renum_impl.F90
)

target_link_libraries(util PRIVATE ${MPI_C_LIBRARIES} ${MPI_Fortran_LIBRARIES} base)
target_link_libraries(util_static PRIVATE ${MPI_C_LIBRARIES} ${MPI_Fortran_LIBRARIES} base)

install(DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  FILES_MATCHING PATTERN "*.mod")

set_target_properties( util_static
  PROPERTIES
  SOVERSION ${PSBLAS_SO_VERSION}
)
set_target_properties ( util
  PROPERTIES
  SOVERSION ${PSBLAS_SO_VERSION}
)
set_target_properties(util_static PROPERTIES OUTPUT_NAME util)
get_target_property(libutil_static util_static OUTPUT_NAME)
set(PSBLAS_LIBS "${CMAKE_INSTALL_LIBDIR}/${CMAKE_STATIC_LIBRARY_PREFIX}${libutil_static}${CMAKE_STATIC_LIBRARY_SUFFIX}")

install(TARGETS util EXPORT PSBLAS_Targets
  DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)

install(TARGETS util_static EXPORT PSBLAS_Targets
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)
