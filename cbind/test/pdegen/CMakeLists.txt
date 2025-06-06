cmake_minimum_required(VERSION 3.10)
project(pdegen C Fortran)

# Check for the installation path for psblas
if(NOT DEFINED PSBLAS_INSTALL_DIR)
    message(FATAL_ERROR "Please specify the path to the psblas installation directory using -DPSBLAS_INSTALL_DIR=<path>")
endif()

# Set the include and library directories based on the provided path
set(INSTALLDIR "${PSBLAS_INSTALL_DIR}")
set(INCDIR "${INSTALLDIR}/include")
set(MODDIR "${INSTALLDIR}/modules")
set(LIBDIR "${INSTALLDIR}/lib")

#set(CMAKE_VERBOSE_MAKEFILE ON) 

# Find the psblas package
find_package(psblas REQUIRED PATHS ${INSTALLDIR})

# Include directories for the Fortran compiler
include_directories(${INCDIR} ${MODDIR})

# Define executable directory
set(EXEDIR "${CMAKE_CURRENT_SOURCE_DIR}/runs")

# Ensure the executable directory exists
file(MAKE_DIRECTORY ${EXEDIR})

# Define sources for the executables
set(SOURCES_PDEGEN pdegen3dc.c)


# Create executables
add_executable(psb_pdegen ${SOURCES_PDEGEN})
target_link_libraries(psb_pdegen psblas::cbind psblas::ext psblas::util psblas::linsolve psblas::prec psblas::base )


# Set output directory for executables
foreach(target psb_pdegen)
    set_target_properties(${target} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${EXEDIR}
    )
endforeach()

