cmake_minimum_required(VERSION 3.10)
project(HelloWorld Fortran)

# Accept a user-defined library path
set(LIBRARY_DIR "" CACHE PATH "Path to the library directory")

# Check if the user provided a library directory
if(NOT LIBRARY_DIR)
    message(FATAL_ERROR "Library directory not specified! Use -DLIBRARY_DIR=path/to/library")
endif()

# Include CMakePackageConfigHelpers to work with configuration files
include(CMakePackageConfigHelpers)

# Find the package
find_package(psblas REQUIRED PATHS ${LIBRARY_DIR}/lib/cmake/psblas NO_DEFAULT_PATH)

# Check if the package was found
if(NOT psblas_FOUND)
    message(FATAL_ERROR "psblas not found!")
endif()

# Include directories for the library
include_directories(${LIBRARY_DIR}/include)  # Path to header files
include_directories(${psblas_DIR}/modules)    # Path to module files

message(STATUS "Library directory: ${psblas_DIR}")

# Add the executables
add_executable(hello hello.f90)
add_executable(pingpong pingpong.f90)

# Link the specific library targets
target_link_libraries(hello PRIVATE psblas::base)
target_link_libraries(pingpong PRIVATE psblas::base)

# Set output directory for executables
set_target_properties(hello pingpong PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/runs)

# Ensure the module path is set correctly
set_target_properties(hello pingpong PROPERTIES
    Fortran_MODULE_DIRECTORY "${psblas_DIR}/modules"
)
