include_directories("${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${mod_dir_tail}")

add_executable( d_matgen
  psb_d_xyz_impl.f90
  psb_d_xyz_mat_mod.f90
  d_matgen.F90 # Main program
)
target_link_libraries(d_matgen base)
