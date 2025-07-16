cmake_minimum_required(VERSION 3.10)
project(psblas_project Fortran)

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

# for Fortran compiler
include_directories(${INCDIR} ${MODDIR})

# Define executable directory
set(EXEDIR "${CMAKE_CURRENT_SOURCE_DIR}/runs")

# Ensure the executable directory exists
file(MAKE_DIRECTORY ${EXEDIR})

# PSBLAS libraries
set(PSBLAS_LIBS  psblas::util psblas::prec psblas::base) 


include(${CMAKE_CURRENT_LIST_DIR}/utils/CMakeLists.txt)

include(${CMAKE_CURRENT_LIST_DIR}/geaxpby/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/gedot/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/spmm/CMakeLists.txt)

# Create executables
add_executable(psb_geaxpby_test ${geaxpby_source_files})
add_executable(psb_gedot_test ${gedot_source_files})
add_executable(psb_spmm_test ${spmm_source_files})

# Link the necessary libraries
target_link_libraries(psb_geaxpby_test ${PSBLAS_LIBS})
target_link_libraries(psb_gedot_test ${PSBLAS_LIBS})
target_link_libraries(psb_spmm_test ${PSBLAS_LIBS})

# Set output directory
set_target_properties(psb_geaxpby_test PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${EXEDIR}
)
set_target_properties(psb_gedot_test PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${EXEDIR}
)
set_target_properties(psb_spmm_test PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${EXEDIR}
)

# Inform the user of successful configuration
message(STATUS "CMake configuration complete. Build with 'make' in the build directory.")
