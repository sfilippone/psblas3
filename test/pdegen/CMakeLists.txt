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

set(EXTLIBDIR "/home/fdurastante/psctoolkit/install/lib")

link_directories(${EXTLIBDIR})

# Include directories for the Fortran compiler
include_directories(${INCDIR} ${MODDIR})

# Define executable directory
set(EXEDIR "${CMAKE_CURRENT_SOURCE_DIR}/runs")

# Ensure the executable directory exists
file(MAKE_DIRECTORY ${EXEDIR})

# Define sources for the executables
set(SOURCES_D_PDE3D psb_d_pde3d.F90)
set(SOURCES_S_PDE3D psb_s_pde3d.F90)
set(SOURCES_D_PDE2D psb_d_pde2d.F90)
set(SOURCES_S_PDE2D psb_s_pde2d.F90)

set(SOURCES_D_PDE_NEST_PSBLAS    psb_d_pde_nest_psblas.F90)
set(SOURCES_D_PDE_NEST_KRYLOV    psb_d_pde_nest_krylov.F90)

# Create executables
add_executable(psb_d_pde3d ${SOURCES_D_PDE3D})
target_link_libraries(psb_d_pde3d psblas::util psblas::linsolve psblas::prec psblas::base)

add_executable(psb_s_pde3d ${SOURCES_S_PDE3D})
target_link_libraries(psb_s_pde3d psblas::util psblas::linsolve psblas::prec psblas::base)

add_executable(psb_d_pde2d ${SOURCES_D_PDE2D})
target_link_libraries(psb_d_pde2d psblas::util psblas::linsolve psblas::prec psblas::base)

add_executable(psb_s_pde2d ${SOURCES_S_PDE2D})
target_link_libraries(psb_s_pde2d psblas::util psblas::linsolve psblas::prec psblas::base)

add_executable(psb_d_pde_nest_psblas ${SOURCES_D_PDE_NEST_PSBLAS})
target_link_libraries(psb_d_pde_nest_psblas psblas::util psblas::linsolve psblas::prec psblas::base)

add_executable(psb_d_pde_nest_krylov ${SOURCES_D_PDE_NEST_KRYLOV})
target_link_libraries(psb_d_pde_nest_krylov psblas::util psblas::linsolve psblas::prec psblas::base)

# Set output directory for executables
foreach(target psb_d_pde3d psb_s_pde3d psb_d_pde2d psb_s_pde2d psb_dist1_pde2d psb_dist2_pde2d psb_dist3_pde2d psb_d_pde_nest_psblas psb_d_pde_nest_krylov)
    set_target_properties(${target} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${EXEDIR}
    )
endforeach()

