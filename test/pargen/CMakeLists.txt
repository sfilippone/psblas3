include_directories("${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${mod_dir_tail}")

add_executable(psb_d_pde2d psb_d_pde2d.f90)
add_executable(psb_d_pde3d psb_d_pde3d.f90)
add_executable(psb_s_pde2d psb_s_pde2d.f90)
add_executable(psb_s_pde3d psb_s_pde3d.f90)

target_link_libraries(psb_d_pde2d krylov util)
target_link_libraries(psb_d_pde3d krylov util)
target_link_libraries(psb_s_pde2d krylov util)
target_link_libraries(psb_s_pde3d krylov util)
