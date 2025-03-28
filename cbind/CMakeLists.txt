set(PSB_cbind_source_files
  base/psb_d_tools_cbind_mod.F90
  base/psb_s_tools_cbind_mod.F90
  base/psb_d_psblas_cbind_mod.f90
  base/psb_objhandle_mod.F90
  base/psb_base_psblas_cbind_mod.f90
  base/psb_z_psblas_cbind_mod.f90
  base/psb_c_comm_cbind_mod.f90
  base/psb_z_serial_cbind_mod.F90
  base/psb_d_serial_cbind_mod.F90
  base/psb_c_tools_cbind_mod.F90
  base/psb_c_serial_cbind_mod.F90
  base/psb_base_string_cbind_mod.f90
  base/psb_base_tools_cbind_mod.F90
  base/psb_z_comm_cbind_mod.f90
  base/psb_s_serial_cbind_mod.F90
  base/psb_base_cbind_mod.f90
  base/psb_s_comm_cbind_mod.f90
  base/psb_s_psblas_cbind_mod.f90
  base/psb_c_psblas_cbind_mod.f90
  base/psb_d_comm_cbind_mod.f90
  base/psb_z_tools_cbind_mod.F90
  base/psb_cpenv_mod.f90
  util/psb_c_util_cbind_mod.f90
  util/psb_s_util_cbind_mod.f90
  util/psb_util_cbind_mod.f90
  util/psb_d_util_cbind_mod.f90
  util/psb_z_util_cbind_mod.f90
  krylov/psb_ckrylov_cbind_mod.f90
  krylov/psb_base_krylov_cbind_mod.f90
  krylov/psb_skrylov_cbind_mod.f90
  krylov/psb_dkrylov_cbind_mod.f90
  krylov/psb_zkrylov_cbind_mod.f90
  prec/psb_dprec_cbind_mod.f90
  prec/psb_cprec_cbind_mod.f90
  prec/psb_prec_cbind_mod.f90
  prec/psb_sprec_cbind_mod.f90
  prec/psb_zprec_cbind_mod.f90
)
foreach(file IN LISTS PSB_cbind_source_files)
  list(APPEND cbind_source_files ${CMAKE_CURRENT_LIST_DIR}/${file})
endforeach()

list(APPEND PSB_cbind_source_C_files
  base/psb_c_dcomm.c
  base/psb_c_scomm.c
  base/psb_c_zcomm.c
  base/psb_c_ccomm.c
  base/psb_c_dbase.c
  base/psb_c_base.c
  base/psb_c_zbase.c
  base/psb_c_cbase.c
  base/psb_c_sbase.c
  prec/psb_c_dprec.c
  prec/psb_c_cprec.c
  prec/psb_c_zprec.c
  prec/psb_c_sprec.c
  #test/pargen/ppdec.c

 )

list(APPEND PSB_cbind_header_C_files
  base/psb_c_sbase.h
  base/psb_c_base.h
  base/psb_c_dcomm.h
  base/psb_c_dbase.h
  base/psb_c_scomm.h
  base/psb_c_ccomm.h
  base/psb_base_cbind.h
  base/psb_c_cbase.h
  base/psb_c_zbase.h
  base/psb_c_zcomm.h
  util/psb_c_zutil.h
  util/psb_c_dutil.h
  util/psb_c_sutil.h
  util/psb_c_cutil.h
  util/psb_util_cbind.h
  krylov/psb_krylov_cbind.h
  prec/psb_c_sprec.h
  prec/psb_c_cprec.h
  prec/psb_prec_cbind.h
  prec/psb_c_dprec.h
  prec/psb_c_zprec.h
)

#if (SERIAL_MPI)
# list(APPEND PSB_base_source_C_files modules/fakempi.c)
#endif()

foreach(file IN LISTS PSB_cbind_source_C_files)
  list(APPEND cbind_source_C_files ${CMAKE_CURRENT_LIST_DIR}/${file})
endforeach()
foreach(file IN LISTS PSB_cbind_header_C_files)
  list(APPEND cbind_header_C_files ${CMAKE_CURRENT_LIST_DIR}/${file})
endforeach()
