include_directories("${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${mod_dir_tail}")

add_executable(psb_cf_sample psb_cf_sample.f90 getp.f90)
add_executable(psb_df_sample psb_df_sample.f90 getp.f90)
add_executable(psb_sf_sample psb_sf_sample.f90 getp.f90)
add_executable(psb_zf_sample psb_zf_sample.f90 getp.f90)

target_link_libraries(psb_cf_sample krylov util base)
target_link_libraries(psb_df_sample krylov util base)
target_link_libraries(psb_sf_sample krylov util base)
target_link_libraries(psb_zf_sample krylov util base)
