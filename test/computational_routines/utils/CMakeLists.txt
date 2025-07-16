set(PSB_gedot_source_files
  psb_test_utils.f90
)
foreach(file IN LISTS PSB_gedot_source_files)
  list(APPEND gedot_source_files ${CMAKE_CURRENT_LIST_DIR}/${file})
endforeach()

