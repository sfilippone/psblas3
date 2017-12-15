include_directories("${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${mod_dir_tail}")

add_executable(hello hello.f90)
add_executable(pingpong pingpong.f90)

target_link_libraries(hello base)
target_link_libraries(pingpong base)
