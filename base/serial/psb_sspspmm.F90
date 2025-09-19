!   
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
!        Salvatore Filippone    
!        Alfredo Buttari      
!   
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
!    
! File:  psb_sspspmm.f90 
! Subroutine: 
! Arguments:
!
!
!
subroutine psb_sspspmm(a,b,c,info)
  use psb_mat_mod
#if !defined(PSB_CMP_INTEL)
  use psb_s_csr_mat_mod
  use psb_s_csc_mat_mod
  use psb_s_serial_mod, psb_protect_name => psb_sspspmm
#endif  
  type(psb_sspmat_type), intent(in)    :: a,b
  type(psb_sspmat_type), intent(out)   :: c
  integer(psb_ipk_), intent(out)                  :: info
  type(psb_s_csr_sparse_mat), allocatable :: ccsr
  type(psb_s_csc_sparse_mat), allocatable :: ccsc
  integer(psb_ipk_) :: err_act
  character(len=*), parameter ::  name='psb_spspmm'
  logical :: done_spmm

#if defined(PSB_CMP_INTEL)  
  interface psb_symbmm
    subroutine psb_ssymbmm(a,b,c,info)
      use psb_s_mat_mod, only : psb_sspmat_type
      import :: psb_ipk_
      implicit none 
      type(psb_sspmat_type), intent(in)  :: a,b
      type(psb_sspmat_type), intent(out) :: c
      integer(psb_ipk_), intent(out)                :: info
    end subroutine psb_ssymbmm
    subroutine psb_sbase_symbmm(a,b,c,info)
      use psb_s_mat_mod, only : psb_s_base_sparse_mat, psb_s_csr_sparse_mat
      import :: psb_ipk_
      implicit none 
      class(psb_s_base_sparse_mat), intent(in) :: a,b
      type(psb_s_csr_sparse_mat), intent(out)  :: c
      integer(psb_ipk_), intent(out)                     :: info
    end subroutine psb_sbase_symbmm
  end interface psb_symbmm

  interface psb_numbmm
    subroutine psb_snumbmm(a,b,c)
      use psb_s_mat_mod, only : psb_sspmat_type
      import :: psb_ipk_
      implicit none 
      type(psb_sspmat_type), intent(in) :: a,b
      type(psb_sspmat_type), intent(inout)  :: c
    end subroutine psb_snumbmm
    subroutine psb_sbase_numbmm(a,b,c)
      use psb_s_mat_mod, only : psb_s_base_sparse_mat, psb_s_csr_sparse_mat
      import :: psb_ipk_
      implicit none 
      class(psb_s_base_sparse_mat), intent(in) :: a,b
      type(psb_s_csr_sparse_mat), intent(inout)  :: c
    end subroutine psb_sbase_numbmm
  end interface psb_numbmm
  interface 
    subroutine psb_scsrspspmm(a,b,c,info)
      use psb_s_mat_mod, only : psb_s_csr_sparse_mat
      import :: psb_ipk_
      implicit none 
      class(psb_s_csr_sparse_mat), intent(in) :: a,b
      type(psb_s_csr_sparse_mat), intent(out) :: c
      integer(psb_ipk_), intent(out)          :: info
    end subroutine psb_scsrspspmm
    subroutine psb_scscspspmm(a,b,c,info)
      use psb_s_mat_mod, only : psb_s_csc_sparse_mat
      import :: psb_ipk_
      implicit none 
      class(psb_s_csc_sparse_mat), intent(in) :: a,b
      type(psb_s_csc_sparse_mat), intent(out) :: c
      integer(psb_ipk_), intent(out)          :: info
    end subroutine psb_scscspspmm
  end interface
#endif
  
  call psb_erractionsave(err_act)
  info = psb_success_

  if ((a%is_null()) .or.(b%is_null())) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif


  ! 
  ! Shortcuts for special cases
  !
  done_spmm = .false. 
  select type(aa=>a%a)
  class is (psb_s_csr_sparse_mat) 
    select type(ba=>b%a)
    class is (psb_s_csr_sparse_mat) 
      
      allocate(ccsr,stat=info)    
      if (info == psb_success_) then 
        call psb_scsrspspmm(aa,ba,ccsr,info)
      else
        info = psb_err_alloc_dealloc_
      end if
      if (info == psb_success_) call move_alloc(ccsr,c%a)
      done_spmm = .true. 

    end select

  class is (psb_s_csc_sparse_mat) 
    select type(ba=>b%a)
    class is (psb_s_csc_sparse_mat) 
      
      allocate(ccsc,stat=info)    
      if (info == psb_success_) then 
        call psb_scscspspmm(aa,ba,ccsc,info)
      else
        info = psb_err_alloc_dealloc_
      end if
      if (info == psb_success_) call move_alloc(ccsc,c%a)
      done_spmm = .true. 

    end select

  end select
  
  !
  ! General code
  !
  if (.not.done_spmm) then 
    call psb_symbmm(a,b,c,info)
    if (info == psb_success_) call psb_numbmm(a,b,c)
  end if
  
  if (info /= psb_success_) then 
    call psb_errpush(info,name) 
    goto 9999
  end if
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_sspspmm

subroutine psb_lsspspmm(a,b,c,info)
  use psb_mat_mod
#if !defined(PSB_CMP_INTEL)
  use psb_s_csr_mat_mod
  use psb_s_csc_mat_mod
  use psb_s_serial_mod, psb_protect_name => psb_lsspspmm
#endif
  implicit none 

  type(psb_lsspmat_type), intent(in)    :: a,b
  type(psb_lsspmat_type), intent(out)   :: c
  integer(psb_ipk_), intent(out)                  :: info
  type(psb_ls_csr_sparse_mat), allocatable :: ccsr
  type(psb_ls_csc_sparse_mat), allocatable :: ccsc
  integer(psb_ipk_) :: err_act
  character(len=*), parameter ::  name='psb_spspmm'
  logical :: done_spmm
#if defined(PSB_CMP_INTEL)
  interface psb_symbmm
    subroutine psb_lssymbmm(a,b,c,info)
      use psb_s_mat_mod, only : psb_lsspmat_type
      import :: psb_ipk_
      implicit none 
      type(psb_lsspmat_type), intent(in)  :: a,b
      type(psb_lsspmat_type), intent(out) :: c
      integer(psb_ipk_), intent(out)                :: info
    end subroutine psb_lssymbmm
    subroutine psb_lsbase_symbmm(a,b,c,info)
      use psb_s_mat_mod, only : psb_ls_base_sparse_mat, psb_ls_csr_sparse_mat
      import :: psb_ipk_
      implicit none 
      class(psb_ls_base_sparse_mat), intent(in) :: a,b
      type(psb_ls_csr_sparse_mat), intent(out)  :: c
      integer(psb_ipk_), intent(out)                     :: info
    end subroutine psb_lsbase_symbmm
  end interface psb_symbmm

  interface psb_numbmm
    subroutine psb_lsnumbmm(a,b,c)
      use psb_s_mat_mod, only : psb_lsspmat_type
      import :: psb_ipk_
      implicit none 
      type(psb_lsspmat_type), intent(in) :: a,b
      type(psb_lsspmat_type), intent(inout)  :: c
    end subroutine psb_lsnumbmm
    subroutine psb_lsbase_numbmm(a,b,c)
      use psb_s_mat_mod, only : psb_ls_base_sparse_mat, psb_ls_csr_sparse_mat
      import :: psb_ipk_
      implicit none 
      class(psb_ls_base_sparse_mat), intent(in) :: a,b
      type(psb_ls_csr_sparse_mat), intent(inout)  :: c
    end subroutine psb_lsbase_numbmm
  end interface psb_numbmm
  interface 
    subroutine psb_lscsrspspmm(a,b,c,info)
      use psb_s_mat_mod, only : psb_ls_csr_sparse_mat
      import :: psb_ipk_
      implicit none 
      class(psb_ls_csr_sparse_mat), intent(in) :: a,b
      type(psb_ls_csr_sparse_mat), intent(out) :: c
      integer(psb_ipk_), intent(out)          :: info
    end subroutine psb_lscsrspspmm
    subroutine psb_lscscspspmm(a,b,c,info)
      use psb_s_mat_mod, only : psb_ls_csc_sparse_mat
      import :: psb_ipk_
      implicit none 
      class(psb_ls_csc_sparse_mat), intent(in) :: a,b
      type(psb_ls_csc_sparse_mat), intent(out) :: c
      integer(psb_ipk_), intent(out)          :: info
    end subroutine psb_lscscspspmm
  end interface
#endif

  call psb_erractionsave(err_act)
  info = psb_success_

  if ((a%is_null()) .or.(b%is_null())) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif


  ! 
  ! Shortcuts for special cases
  !
  done_spmm = .false. 
  select type(aa=>a%a)
  class is (psb_ls_csr_sparse_mat) 
    select type(ba=>b%a)
    class is (psb_ls_csr_sparse_mat) 
      
      allocate(ccsr,stat=info)    
      if (info == psb_success_) then 
        call psb_lscsrspspmm(aa,ba,ccsr,info)
      else
        info = psb_err_alloc_dealloc_
      end if
      if (info == psb_success_) call move_alloc(ccsr,c%a)
      done_spmm = .true. 

    end select

  class is (psb_ls_csc_sparse_mat) 
    select type(ba=>b%a)
    class is (psb_ls_csc_sparse_mat) 
      
      allocate(ccsc,stat=info)    
      if (info == psb_success_) then 
        call psb_lscscspspmm(aa,ba,ccsc,info)
      else
        info = psb_err_alloc_dealloc_
      end if
      if (info == psb_success_) call move_alloc(ccsc,c%a)
      done_spmm = .true. 

    end select

  end select
  
  !
  ! General code
  !
  if (.not.done_spmm) then 
    call psb_symbmm(a,b,c,info)
    if (info == psb_success_) call psb_numbmm(a,b,c)
  end if
  
  if (info /= psb_success_) then 
    call psb_errpush(info,name) 
    goto 9999
  end if
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_lsspspmm

