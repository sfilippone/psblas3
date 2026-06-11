cmake_minimum_required(VERSION 3.10)
project(nested Fortran)

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
file(MAKE_DIRECTORY ${EXEDIR})

# Nested (block-structured / MATNEST) tests
set(SOURCES_D_NEST_GLOB_TEST  psb_d_nest_glob_test.F90)
set(SOURCES_D_NEST_RECT_TEST  psb_d_nest_rect_test.F90)
set(SOURCES_D_NEST_CG_TEST    psb_d_nest_cg_test.F90)
set(SOURCES_D_NEST_BUILDER_TEST psb_d_nest_builder_test.F90)

add_executable(psb_d_nest_glob_test ${SOURCES_D_NEST_GLOB_TEST})
target_link_libraries(psb_d_nest_glob_test psblas::util psblas::linsolve psblas::prec psblas::base)

add_executable(psb_d_nest_rect_test ${SOURCES_D_NEST_RECT_TEST})
target_link_libraries(psb_d_nest_rect_test psblas::util psblas::linsolve psblas::prec psblas::base)

add_executable(psb_d_nest_cg_test ${SOURCES_D_NEST_CG_TEST})
target_link_libraries(psb_d_nest_cg_test psblas::util psblas::linsolve psblas::prec psblas::base)

add_executable(psb_d_nest_builder_test ${SOURCES_D_NEST_BUILDER_TEST})
target_link_libraries(psb_d_nest_builder_test psblas::util psblas::linsolve psblas::prec psblas::base)

# Set output directory for executables
foreach(target psb_d_nest_glob_test psb_d_nest_rect_test psb_d_nest_cg_test psb_d_nest_builder_test)
    set_target_properties(${target} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${EXEDIR}
    )
endforeach()
