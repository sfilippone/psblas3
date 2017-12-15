include_directories("${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${mod_dir_tail}")

add_executable(fortran_interface_test
  psb_c_mvsv_tester.f90
  psb_d_mvsv_tester.f90
  psb_mvsv_tester.f90
  psb_s_mvsv_tester.f90
  psbtf.f90 # Main program
  psb_z_mvsv_tester.f90
)
target_link_libraries(fortran_interface_test base)
