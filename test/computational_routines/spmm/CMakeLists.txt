set(PSB_spmm_source_files
  psb_spmm_test.f90
  spmm.f90
)
foreach(file IN LISTS PSB_spmm_source_files)
  list(APPEND spmm_source_files ${CMAKE_CURRENT_LIST_DIR}/${file})
endforeach()

