include_directories("${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${mod_dir_tail}")

add_executable(d_file_spmv  d_file_spmv.f90)
add_executable(pdgenspmv    pdgenspmv.f90)
add_executable(s_file_spmv  s_file_spmv.f90)

target_link_libraries(d_file_spmv  util)
target_link_libraries(pdgenspmv    util)
target_link_libraries(s_file_spmv  util)
