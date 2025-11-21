cmake_minimum_required(VERSION 3.10)
project(pargen Fortran)

# Check for the installation path for psblas
if(NOT DEFINED PSBLAS_INSTALL_DIR)
    message(FATAL_ERROR "Please specify the path to the psblas installation directory using -DPSBLAS_INSTALL_DIR=<path>")
endif()

# Set the include and library directories based on the provided path
set(INSTALLDIR "${PSBLAS_INSTALL_DIR}")
set(INCDIR "${INSTALLDIR}/include")
set(MODDIR "${INSTALLDIR}/modules")
set(LIBDIR "${INSTALLDIR}/lib")

# Find the psblas package
find_package(psblas REQUIRED PATHS ${INSTALLDIR})

# Include directories for the Fortran compiler
include_directories(${INCDIR} ${MODDIR})

# Define executable directory
set(EXEDIR "${CMAKE_CURRENT_SOURCE_DIR}/runs")

# Ensure the executable directory exists
file(MAKE_DIRECTORY ${EXEDIR})

# Source file use in all the targets
set(COMMON_SOURCE getp.f90)


# Define sources for the executables
set(SOURCES_SFOBJS ${COMMON_SOURCE} psb_sf_sample.f90)
set(SOURCES_DFOBJS ${COMMON_SOURCE} psb_df_sample.f90)
set(SOURCES_CFOBJS ${COMMON_SOURCE} psb_cf_sample.f90)
set(SOURCES_ZFOBJS ${COMMON_SOURCE} psb_zf_sample.f90)

# Create executables
add_executable(psb_sf_sample ${SOURCES_SFOBJS})
target_link_libraries(psb_sf_sample psblas::util psblas::linsolve psblas::prec psblas::base)

add_executable(psb_df_sample ${SOURCES_DFOBJS})
target_link_libraries(psb_df_sample psblas::util psblas::linsolve psblas::prec psblas::base)

add_executable(psb_cf_sample ${SOURCES_CFOBJS})
target_link_libraries(psb_cf_sample psblas::util psblas::linsolve psblas::prec psblas::base)

add_executable(psb_zf_sample ${SOURCES_ZFOBJS})
target_link_libraries(psb_zf_sample psblas::util psblas::linsolve psblas::prec psblas::base)

# Set output directory for executables
foreach(target psb_sf_sample psb_df_sample psb_cf_sample psb_zf_sample)
    set_target_properties(${target} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${EXEDIR}
    )
endforeach()

