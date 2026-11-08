cmake_minimum_required(VERSION 3.16)
project(psblas_c_sparse_ops C Fortran)

find_package(MPI REQUIRED COMPONENTS C Fortran)
find_package(psblas CONFIG REQUIRED HINTS "${PSBLAS_INSTALL_DIR}")
enable_testing()

foreach(family d s c z)
  set(target ${family}_sparse_ops)
  add_executable(${target} sparse_ops.c)
  set_target_properties(${target} PROPERTIES
    C_STANDARD 99 C_STANDARD_REQUIRED YES LINKER_LANGUAGE Fortran)
  target_link_libraries(${target} PRIVATE
    psblas::cbind psblas::util psblas::linsolve psblas::prec
    psblas::ext psblas::base MPI::MPI_C MPI::MPI_Fortran m)
  if(family STREQUAL "s")
    target_compile_definitions(${target} PRIVATE TEST_SINGLE_REAL)
  elseif(family STREQUAL "c")
    target_compile_definitions(${target} PRIVATE TEST_SINGLE_COMPLEX)
  elseif(family STREQUAL "z")
    target_compile_definitions(${target} PRIVATE TEST_DOUBLE_COMPLEX)
  endif()
  foreach(ranks 1 2 4)
    add_test(NAME ${target}_${ranks}_ranks COMMAND
      ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks}
      ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${target}> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(${target}_${ranks}_ranks PROPERTIES
      PROCESSORS ${ranks} TIMEOUT 60)
  endforeach()
endforeach()
