!
! Parallel Sparse BLAS fortran interface testing code
!
! 
!

program main

  use psb_base_mod
  use psb_mvsv_tester
  implicit none
  integer(psb_ipk_), parameter :: psb_fidasize_=16
  integer(psb_ipk_) :: res,passed=0,failed=0;
  integer(psb_ipk_) :: ictxt, iam=-1, np=-1
  character(len=psb_fidasize_) :: afmt

  write(psb_out_unit,*) 'Format ?'
  read(psb_inp_unit,*) afmt
!  afmt = 'COO'

  call psb_init(ictxt)
  call psb_info(ictxt,iam,np)
  if(iam<0)then
    goto 9999
  endif
  call       s_usmv_2_n_ap3_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_usmv_2_t_ap3_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_usmv_2_c_ap3_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_usmv_2_n_ap3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_usmv_2_t_ap3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_usmv_2_c_ap3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_usmv_2_n_ap1_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_usmv_2_t_ap1_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_usmv_2_c_ap1_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_usmv_2_n_ap1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_usmv_2_t_ap1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_usmv_2_c_ap1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_usmv_2_n_am1_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_usmv_2_t_am1_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_usmv_2_c_am1_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_usmv_2_n_am1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_usmv_2_t_am1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_usmv_2_c_am1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_usmv_2_n_am3_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_usmv_2_t_am3_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_usmv_2_c_am3_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_usmv_2_n_am3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_usmv_2_t_am3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_usmv_2_c_am3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_ussv_2_n_ap3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_ussv_2_t_ap3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_ussv_2_c_ap3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_ussv_2_n_ap1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_ussv_2_t_ap1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_ussv_2_c_ap1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_ussv_2_n_am1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_ussv_2_t_am1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_ussv_2_c_am1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_ussv_2_n_am3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_ussv_2_t_am3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       s_ussv_2_c_am3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_usmv_2_n_ap3_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_usmv_2_t_ap3_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_usmv_2_c_ap3_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_usmv_2_n_ap3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_usmv_2_t_ap3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_usmv_2_c_ap3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_usmv_2_n_ap1_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_usmv_2_t_ap1_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_usmv_2_c_ap1_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_usmv_2_n_ap1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_usmv_2_t_ap1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_usmv_2_c_ap1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_usmv_2_n_am1_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_usmv_2_t_am1_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_usmv_2_c_am1_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_usmv_2_n_am1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_usmv_2_t_am1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_usmv_2_c_am1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_usmv_2_n_am3_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_usmv_2_t_am3_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_usmv_2_c_am3_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_usmv_2_n_am3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_usmv_2_t_am3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_usmv_2_c_am3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_ussv_2_n_ap3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_ussv_2_t_ap3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_ussv_2_c_ap3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_ussv_2_n_ap1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_ussv_2_t_ap1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_ussv_2_c_ap1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_ussv_2_n_am1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_ussv_2_t_am1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_ussv_2_c_am1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_ussv_2_n_am3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_ussv_2_t_am3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       d_ussv_2_c_am3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_usmv_2_n_ap3_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_usmv_2_t_ap3_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_usmv_2_c_ap3_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_usmv_2_n_ap3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_usmv_2_t_ap3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_usmv_2_c_ap3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_usmv_2_n_ap1_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_usmv_2_t_ap1_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_usmv_2_c_ap1_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_usmv_2_n_ap1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_usmv_2_t_ap1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_usmv_2_c_ap1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_usmv_2_n_am1_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_usmv_2_t_am1_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_usmv_2_c_am1_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_usmv_2_n_am1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_usmv_2_t_am1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_usmv_2_c_am1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_usmv_2_n_am3_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_usmv_2_t_am3_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_usmv_2_c_am3_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_usmv_2_n_am3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_usmv_2_t_am3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_usmv_2_c_am3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_ussv_2_n_ap3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_ussv_2_t_ap3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_ussv_2_c_ap3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_ussv_2_n_ap1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_ussv_2_t_ap1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_ussv_2_c_ap1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_ussv_2_n_am1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_ussv_2_t_am1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_ussv_2_c_am1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_ussv_2_n_am3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_ussv_2_t_am3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       c_ussv_2_c_am3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_usmv_2_n_ap3_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_usmv_2_t_ap3_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_usmv_2_c_ap3_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_usmv_2_n_ap3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_usmv_2_t_ap3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_usmv_2_c_ap3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_usmv_2_n_ap1_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_usmv_2_t_ap1_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_usmv_2_c_ap1_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_usmv_2_n_ap1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_usmv_2_t_ap1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_usmv_2_c_ap1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_usmv_2_n_am1_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_usmv_2_t_am1_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_usmv_2_c_am1_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_usmv_2_n_am1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_usmv_2_t_am1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_usmv_2_c_am1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_usmv_2_n_am3_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_usmv_2_t_am3_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_usmv_2_c_am3_bp1_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_usmv_2_n_am3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_usmv_2_t_am3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_usmv_2_c_am3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_ussv_2_n_ap3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_ussv_2_t_ap3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_ussv_2_c_ap3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_ussv_2_n_ap1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_ussv_2_t_ap1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_ussv_2_c_ap1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_ussv_2_n_am1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_ussv_2_t_am1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_ussv_2_c_am1_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_ussv_2_n_am3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_ussv_2_t_am3_bm0_ix1_iy1(res,afmt,ictxt)
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

  call       z_ussv_2_c_am3_bm0_ix1_iy1(res,afmt,ictxt)        
  if(res /= 0)failed=failed+1
  if(res.eq.0)passed=passed+1
  res=0

9999 continue
  print *,"PASSED:",passed
  print *,"FAILED:",failed
  call psb_exit(ictxt)

end program main



