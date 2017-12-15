include_directories("${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${mod_dir_tail}")

add_executable(dhb2mm dhb2mm.f90)
add_executable(dmm2hb dmm2hb.f90)
add_executable(zhb2mm zhb2mm.f90)
add_executable(zmm2hb zmm2hb.f90)

target_link_libraries(dhb2mm util)
target_link_libraries(dmm2hb util)
target_link_libraries(zhb2mm util)
target_link_libraries(zmm2hb util)
