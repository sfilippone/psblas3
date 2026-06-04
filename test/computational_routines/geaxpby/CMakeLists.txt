set(PSB_geaxpby_source_files
  psb_geaxpby_test.f90
  geaxpby.f90
)
foreach(file IN LISTS PSB_geaxpby_source_files)
  list(APPEND geaxpby_source_files ${CMAKE_CURRENT_LIST_DIR}/${file})
endforeach()

