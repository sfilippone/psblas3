module psb_s_bjacprec

  use psb_s_base_prec_mod
  
  type, extends(psb_s_base_prec_type)   :: psb_s_bjac_prec_type
    integer, allocatable                :: iprcparm(:)
    type(psb_sspmat_type), allocatable  :: av(:)
    type(psb_s_vect_type), allocatable  :: dv
  contains
    procedure, pass(prec) :: s_apply_v => psb_s_bjac_apply_vect
    procedure, pass(prec) :: s_apply   => psb_s_bjac_apply
    procedure, pass(prec) :: precbld   => psb_s_bjac_precbld
    procedure, pass(prec) :: precinit  => psb_s_bjac_precinit
    procedure, pass(prec) :: precseti  => psb_s_bjac_precseti
    procedure, pass(prec) :: precsetr  => psb_s_bjac_precsetr
    procedure, pass(prec) :: precsetc  => psb_s_bjac_precsetc
    procedure, pass(prec) :: precfree  => psb_s_bjac_precfree
    procedure, pass(prec) :: precdescr => psb_s_bjac_precdescr
    procedure, pass(prec) :: dump      => psb_s_bjac_dump
    procedure, pass(prec) :: sizeof    => psb_s_bjac_sizeof
    procedure, pass(prec) :: get_nzeros => psb_s_bjac_get_nzeros
  end type psb_s_bjac_prec_type

  private :: psb_s_bjac_sizeof, psb_s_bjac_precdescr, psb_s_bjac_get_nzeros
 

  character(len=15), parameter, private :: &
       &  fact_names(0:2)=(/'None          ','ILU(n)        ',&
       &  'ILU(eps)      '/)

  
  interface  
    subroutine psb_s_bjac_dump(prec,info,prefix,head)
      import :: psb_desc_type, psb_s_bjac_prec_type, psb_s_vect_type, psb_spk_
      class(psb_s_bjac_prec_type), intent(in) :: prec
      integer, intent(out)                    :: info
      character(len=*), intent(in), optional  :: prefix,head
    end subroutine psb_s_bjac_dump
  end interface

  interface  
    subroutine psb_s_bjac_apply_vect(alpha,prec,x,beta,y,desc_data,info,trans,work)
      import :: psb_desc_type, psb_s_bjac_prec_type, psb_s_vect_type, psb_spk_
      type(psb_desc_type),intent(in)    :: desc_data
      class(psb_s_bjac_prec_type), intent(inout)  :: prec
      real(psb_spk_),intent(in)         :: alpha,beta
      type(psb_s_vect_type),intent(inout)   :: x
      type(psb_s_vect_type),intent(inout)   :: y
      integer, intent(out)              :: info
      character(len=1), optional        :: trans
      real(psb_spk_),intent(inout), optional, target :: work(:)
    end subroutine psb_s_bjac_apply_vect
  end interface

  interface
    subroutine psb_s_bjac_apply(alpha,prec,x,beta,y,desc_data,info,trans,work)
      import :: psb_desc_type, psb_s_bjac_prec_type, psb_s_vect_type, psb_spk_
      
      type(psb_desc_type),intent(in)    :: desc_data
      class(psb_s_bjac_prec_type), intent(in)  :: prec
      real(psb_spk_),intent(in)         :: alpha,beta
      real(psb_spk_),intent(inout)      :: x(:)
      real(psb_spk_),intent(inout)      :: y(:)
      integer, intent(out)              :: info
      character(len=1), optional        :: trans
      real(psb_spk_),intent(inout), optional, target :: work(:)
    end subroutine psb_s_bjac_apply
  end interface
  
  interface
    subroutine psb_s_bjac_precinit(prec,info)
      import :: psb_desc_type, psb_s_bjac_prec_type, psb_s_vect_type, psb_spk_
      class(psb_s_bjac_prec_type),intent(inout) :: prec
      integer, intent(out)                     :: info
    end subroutine psb_s_bjac_precinit
  end interface
  
  interface
    subroutine psb_s_bjac_precbld(a,desc_a,prec,info,upd,amold,afmt,vmold)
      import :: psb_desc_type, psb_s_bjac_prec_type, psb_s_vect_type, psb_spk_, &
           & psb_sspmat_type, psb_s_base_sparse_mat, psb_s_base_vect_type
      type(psb_sspmat_type), intent(in), target :: a
      type(psb_desc_type), intent(in), target   :: desc_a
      class(psb_s_bjac_prec_type),intent(inout) :: prec
      integer, intent(out)                      :: info
      character, intent(in), optional           :: upd
      character(len=*), intent(in), optional    :: afmt
      class(psb_s_base_sparse_mat), intent(in), optional :: amold
      class(psb_s_base_vect_type), intent(in), optional  :: vmold
    end subroutine psb_s_bjac_precbld
  end interface
  
  interface
    subroutine psb_s_bjac_precseti(prec,what,val,info)
      import :: psb_desc_type, psb_s_bjac_prec_type, psb_s_vect_type, psb_spk_
      class(psb_s_bjac_prec_type),intent(inout) :: prec
      integer, intent(in)                      :: what 
      integer, intent(in)                      :: val 
      integer, intent(out)                     :: info
    end subroutine psb_s_bjac_precseti
  end interface
  
!!$  interface
!!$    subroutine psb_s_bjac_precsetr(prec,what,val,info)
!!$      import :: psb_desc_type, psb_s_bjac_prec_type, psb_s_vect_type, psb_spk_
!!$      class(psb_s_bjac_prec_type),intent(inout) :: prec
!!$      integer, intent(in)                      :: what 
!!$      real(psb_spk_), intent(in)               :: val 
!!$      integer, intent(out)                     :: info
!!$    end subroutine psb_s_bjac_precsetr
!!$  end interface
!!$  
!!$  interface
!!$    subroutine psb_s_bjac_precsetc(prec,what,val,info)
!!$      import :: psb_desc_type, psb_s_bjac_prec_type, psb_s_vect_type, psb_spk_
!!$      class(psb_s_bjac_prec_type),intent(inout) :: prec
!!$      integer, intent(in)                      :: what 
!!$      character(len=*), intent(in)             :: val
!!$      integer, intent(out)                     :: info
!!$    end subroutine psb_s_bjac_precsetc
!!$  end interface
!!$  
!!$  interface 
!!$    subroutine psb_s_bjac_precfree(prec,info)
!!$      import :: psb_desc_type, psb_s_bjac_prec_type, psb_s_vect_type, psb_spk_
!!$      class(psb_s_bjac_prec_type), intent(inout) :: prec
!!$      integer, intent(out)                :: info
!!$    end subroutine psb_s_bjac_precfree
!!$  end interface 

contains

  subroutine psb_s_bjac_precdescr(prec,iout)
    
    Implicit None

    class(psb_s_bjac_prec_type), intent(in) :: prec
    integer, intent(in), optional    :: iout

    Integer :: err_act, nrow, info
    character(len=20)  :: name='s_bjac_precdescr'
    integer :: iout_

    call psb_erractionsave(err_act)

    info = psb_success_
   
    if (present(iout)) then 
      iout_ = iout
    else
      iout_ = 6 
    end if

    if (.not.allocated(prec%iprcparm)) then 
      info = 1124
      call psb_errpush(info,name,a_err="preconditioner")
      goto 9999
    end if
    
    write(iout_,*) 'Block Jacobi with: ',&
         &  fact_names(prec%iprcparm(psb_f_type_))
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine psb_s_bjac_precdescr


  function psb_s_bjac_sizeof(prec) result(val)
    class(psb_s_bjac_prec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    
    val = 0
    if (allocated(prec%dv)) then 
      val = val + psb_sizeof_sp * prec%dv%get_nrows()
    endif
    if (allocated(prec%av)) then 
      val = val + prec%av(psb_l_pr_)%sizeof()
      val = val + prec%av(psb_u_pr_)%sizeof()
    endif
    return
  end function psb_s_bjac_sizeof

  function psb_s_bjac_get_nzeros(prec) result(val)

    class(psb_s_bjac_prec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    
    val = 0
    if (allocated(prec%dv)) then 
      val = val + prec%dv%get_nrows()
    endif
    if (allocated(prec%av)) then 
      val = val + prec%av(psb_l_pr_)%get_nzeros()
      val = val + prec%av(psb_u_pr_)%get_nzeros()
    endif
    return
  end function psb_s_bjac_get_nzeros


  subroutine psb_s_bjac_precsetr(prec,what,val,info)

    Implicit None

    class(psb_s_bjac_prec_type),intent(inout) :: prec
    integer, intent(in)                      :: what 
    real(psb_spk_), intent(in)               :: val 
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='s_bjac_precset'

    call psb_erractionsave(err_act)

    info = psb_success_

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine psb_s_bjac_precsetr

  subroutine psb_s_bjac_precsetc(prec,what,val,info)

    Implicit None

    class(psb_s_bjac_prec_type),intent(inout) :: prec
    integer, intent(in)                      :: what 
    character(len=*), intent(in)             :: val
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='s_bjac_precset'

    call psb_erractionsave(err_act)

    info = psb_success_

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine psb_s_bjac_precsetc

  subroutine psb_s_bjac_precfree(prec,info)

    Implicit None

    class(psb_s_bjac_prec_type), intent(inout) :: prec
    integer, intent(out)                :: info

    Integer :: err_act, i
    character(len=20)  :: name='s_bjac_precfree'

    call psb_erractionsave(err_act)

    info = psb_success_
    if (allocated(prec%av)) then 
      do i=1,size(prec%av) 
        call prec%av(i)%free()
      enddo
      deallocate(prec%av,stat=info)
    end if

    if (allocated(prec%dv)) then 
      call prec%dv%free(info)
      if (info == 0) deallocate(prec%dv,stat=info)
    end if
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine psb_s_bjac_precfree

end module psb_s_bjacprec
