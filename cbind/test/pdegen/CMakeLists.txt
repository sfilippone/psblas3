cmake_minimum_required(VERSION 3.15)
project(PSB_PDEGen3D C)

# Project structure (mimicking Makefile's TOP)
set(TOP_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../..") # Adjust if needed

# Include directories
include_directories(
    "."
    "${CMAKE_CURRENT_SOURCE_DIR}/../.."
    "${TOP_DIR}/include"
)

# Fortran module directory
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")
set(FMFLAG "-I") # Adjust if needed

# Library directories
link_directories("${TOP_DIR}/lib")

# Libraries
set(PSBC_LIBS psb_cbind)
set(PSB_LIBS psb_util psb_linsolve psb_prec psb_base)

# Executable directory
set(EXEDIR "${CMAKE_CURRENT_BINARY_DIR}/runs")
file(MAKE_DIRECTORY "${EXEDIR}")

# Source files
set(SOURCE_FILES pdegen3dc.c)

# Create executable
add_executable(pdegen3dc ${SOURCE_FILES})

# Link libraries
target_link_libraries(pdegen3dc PUBLIC
    ${PSBC_LIBS}
    ${PSB_LIBS}
    m  # Math library
    gfortran # gfortran library
)

# Move executable to EXEDIR (post-build)
add_custom_command(TARGET pdegen3dc POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E move $<TARGET_FILE:pdegen3dc> ${EXEDIR}
)

# Create "runs" directory (if it doesn't exist)
add_custom_target(create_runs_dir
    COMMAND ${CMAKE_COMMAND} -E make_directory "${EXEDIR}"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Creating runs directory"
)
add_dependencies(pdegen3dc create_runs_dir)

# verycleanlib target
add_custom_target(verycleanlib
    COMMAND ${CMAKE_COMMAND} -E chdir "../../" ${CMAKE_COMMAND} --build . --target veryclean
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Cleaning library"
)

# lib target
add_custom_target(lib
    COMMAND ${CMAKE_COMMAND} -E chdir "../../" ${CMAKE_COMMAND} --build . --target library
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Building library"
)

# tests target
add_custom_target(tests
    COMMAND ${CMAKE_COMMAND} -E chdir "${EXEDIR}" ${CMAKE_CURRENT_SOURCE_DIR}/${EXEDIR}/pdegen3dc < ${CMAKE_CURRENT_SOURCE_DIR}/${EXEDIR}/pdegen3d.inp
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS pdegen3dc
    COMMENT "Running tests"
)

add_dependencies(tests pdegen3dc)

# Install target (optional)
install(TARGETS pdegen3dc DESTINATION bin)
