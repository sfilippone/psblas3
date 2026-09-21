cmake_minimum_required(VERSION 3.15)

# This directory can also be built independently against an installed PSBLAS.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(PSBLASNestedTests Fortran)
  find_package(psblas REQUIRED HINTS ${PSBLAS_INSTALL_DIR})
  find_package(MPI REQUIRED COMPONENTS Fortran)
  set(NEST_LIBS psblas::util psblas::linsolve psblas::prec psblas::ext psblas::base MPI::MPI_Fortran)
  if(PSBLAS_INSTALL_DIR)
    include_directories("${PSBLAS_INSTALL_DIR}/include" "${PSBLAS_INSTALL_DIR}/modules")
  endif()
else()
  set(NEST_LIBS util linsolve prec ext base MPI::MPI_Fortran)
  include_directories("${PROJECT_BINARY_DIR}/include" "${PROJECT_BINARY_DIR}/modules")
endif()

# Executable directory
set(EXEDIR "${CMAKE_CURRENT_BINARY_DIR}/runs")

# Nested (block-structured / MATNEST) tests
set(SOURCES_D_NEST_GLOB_TEST psb_d_nest_glob_test.F90)
set(SOURCES_D_NEST_RECT_TEST psb_d_nest_rect_test.F90)
set(SOURCES_D_NEST_CG_TEST psb_d_nest_cg_test.F90)
set(SOURCES_D_NEST_PREC_CG_TEST psb_d_nest_prec_cg_test.F90)
set(SOURCES_D_NEST_MULT_PREC_CG_TEST psb_d_nest_mult_prec_cg_test.F90)
set(SOURCES_D_NEST_SCHUR_PREC_CG_TEST psb_d_nest_schur_prec_cg_test.F90)
set(SOURCES_D_NEST_PETSC_STOKES_TEST psb_d_nest_petsc_stokes_test.F90)

add_executable(psb_d_nest_glob_test ${SOURCES_D_NEST_GLOB_TEST})
target_link_libraries(psb_d_nest_glob_test PRIVATE ${NEST_LIBS})
set_target_properties(psb_d_nest_glob_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${EXEDIR}")

add_executable(psb_d_nest_rect_test ${SOURCES_D_NEST_RECT_TEST})
target_link_libraries(psb_d_nest_rect_test PRIVATE ${NEST_LIBS})
set_target_properties(psb_d_nest_rect_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${EXEDIR}")

add_executable(psb_d_nest_cg_test ${SOURCES_D_NEST_CG_TEST})
target_link_libraries(psb_d_nest_cg_test PRIVATE ${NEST_LIBS})
set_target_properties(psb_d_nest_cg_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${EXEDIR}")

add_executable(psb_d_nest_prec_cg_test ${SOURCES_D_NEST_PREC_CG_TEST})
target_link_libraries(psb_d_nest_prec_cg_test PRIVATE ${NEST_LIBS})
set_target_properties(psb_d_nest_prec_cg_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${EXEDIR}")

add_executable(psb_d_nest_mult_prec_cg_test ${SOURCES_D_NEST_MULT_PREC_CG_TEST})
target_link_libraries(psb_d_nest_mult_prec_cg_test PRIVATE ${NEST_LIBS})
set_target_properties(psb_d_nest_mult_prec_cg_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${EXEDIR}")

add_executable(psb_d_nest_schur_prec_cg_test ${SOURCES_D_NEST_SCHUR_PREC_CG_TEST})
target_link_libraries(psb_d_nest_schur_prec_cg_test PRIVATE ${NEST_LIBS})
set_target_properties(psb_d_nest_schur_prec_cg_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${EXEDIR}")

add_executable(psb_d_nest_petsc_stokes_test ${SOURCES_D_NEST_PETSC_STOKES_TEST})
target_link_libraries(psb_d_nest_petsc_stokes_test PRIVATE ${NEST_LIBS})
set_target_properties(psb_d_nest_petsc_stokes_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${EXEDIR}")