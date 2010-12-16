module psb_indx_map_mod
  use psb_const_mod
  use psb_desc_const_mod
  
  type      :: psb_indx_map

    integer :: state          = psb_desc_null_
    integer :: ictxt          = -1
    integer :: mpic           = -1
    integer :: global_rows    = -1
    integer :: global_cols    = -1
    integer :: local_rows     = -1
    integer :: local_cols     = -1

  contains

    procedure, pass(idxmap)  :: get_state => base_get_state
    procedure, pass(idxmap)  :: set_state => base_set_state
    procedure, pass(idxmap)  :: is_null   => base_is_null
    procedure, pass(idxmap)  :: is_repl   => base_is_repl
    procedure, pass(idxmap)  :: is_bld    => base_is_bld
    procedure, pass(idxmap)  :: is_upd    => base_is_upd
    procedure, pass(idxmap)  :: is_asb    => base_is_asb
    procedure, pass(idxmap)  :: is_valid  => base_is_valid
    procedure, pass(idxmap)  :: is_ovl    => base_is_ovl
    procedure, pass(idxmap)  :: get_gr    => base_get_gr
    procedure, pass(idxmap)  :: get_gc    => base_get_gc
    procedure, pass(idxmap)  :: get_lr    => base_get_lr
    procedure, pass(idxmap)  :: get_lc    => base_get_lc
    procedure, pass(idxmap)  :: get_ctxt  => base_get_ctxt
    procedure, pass(idxmap)  :: get_mpic  => base_get_mpic
    procedure, pass(idxmap)  :: sizeof    => base_sizeof
    procedure, pass(idxmap)  :: set_null  => base_set_null
    procedure, pass(idxmap)  :: row_extendable => base_row_extendable

    procedure, pass(idxmap)  :: set_gr    => base_set_gr
    procedure, pass(idxmap)  :: set_gc    => base_set_gc
    procedure, pass(idxmap)  :: set_lr    => base_set_lr
    procedure, pass(idxmap)  :: set_lc    => base_set_lc
    procedure, pass(idxmap)  :: set_ctxt  => base_set_ctxt
    procedure, pass(idxmap)  :: set_mpic  => base_set_mpic
    
    procedure, pass(idxmap)  :: get_fmt   => base_get_fmt

    procedure, pass(idxmap)  :: asb   => base_asb
    procedure, pass(idxmap)  :: free  => base_free

    procedure, pass(idxmap)  :: l2gs1  => base_l2gs1
    procedure, pass(idxmap)  :: l2gs2  => base_l2gs2
    procedure, pass(idxmap)  :: l2gv1  => base_l2gv1
    procedure, pass(idxmap)  :: l2gv2  => base_l2gv2
    generic, public          :: l2g => l2gs1, l2gs2, l2gv1, l2gv2

    procedure, pass(idxmap)  :: g2ls1  => base_g2ls1
    procedure, pass(idxmap)  :: g2ls2  => base_g2ls2
    procedure, pass(idxmap)  :: g2lv1  => base_g2lv1
    procedure, pass(idxmap)  :: g2lv2  => base_g2lv2
    generic, public          :: g2l => g2ls1, g2ls2, g2lv1, g2lv2

    procedure, pass(idxmap)  :: g2ls1_ins  => base_g2ls1_ins
    procedure, pass(idxmap)  :: g2ls2_ins  => base_g2ls2_ins
    procedure, pass(idxmap)  :: g2lv1_ins  => base_g2lv1_ins
    procedure, pass(idxmap)  :: g2lv2_ins  => base_g2lv2_ins
    generic, public          :: g2l_ins => g2ls1_ins, g2ls2_ins,&
         &                     g2lv1_ins, g2lv2_ins

    procedure, pass(idxmap)  :: fnd_owner => psb_indx_map_fnd_owner
    procedure, pass(idxmap)  :: init_vl   => base_init_vl
    generic, public          :: init      => init_vl

  end type psb_indx_map

  private :: base_get_state, base_set_state, base_is_repl, base_is_bld,&
       & base_is_upd, base_is_asb, base_is_valid, base_is_ovl,&
       & base_get_gr, base_get_gc, base_get_lr, base_get_lc, base_get_ctxt,&
       & base_get_mpic, base_sizeof, base_set_null, base_set_gr,&
       & base_set_gc, base_set_lr, base_set_lc, base_set_ctxt,&
       & base_set_mpic, base_get_fmt, base_asb, base_free,&
       & base_l2gs1, base_l2gs2, base_l2gv1, base_l2gv2,&
       & base_g2ls1, base_g2ls2, base_g2lv1, base_g2lv2,&
       & base_g2ls1_ins, base_g2ls2_ins, base_g2lv1_ins,&
       & base_g2lv2_ins, base_init_vl, base_is_null, base_row_extendable

  interface 
    subroutine psb_indx_map_fnd_owner(idx,iprc,idxmap,info)
      import :: psb_indx_map
      implicit none 
      integer, intent(in) :: idx(:)
      integer, allocatable, intent(out) ::  iprc(:)
      class(psb_indx_map), intent(in) :: idxmap
      integer, intent(out) :: info
    end subroutine psb_indx_map_fnd_owner
  end interface
  
contains

  function base_get_state(idxmap) result(val)
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer :: val
    
    val = idxmap%state

  end function base_get_state
  

  function base_get_gr(idxmap) result(val)
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer :: val
    
    val = idxmap%global_rows

  end function base_get_gr
  

  function base_get_gc(idxmap) result(val)
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer :: val
    
    val = idxmap%global_cols

  end function base_get_gc
  

  function base_get_lr(idxmap) result(val)
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer :: val
    
    val = idxmap%local_rows

  end function base_get_lr
  

  function base_get_lc(idxmap) result(val)
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer :: val
    
    val = idxmap%local_cols

  end function base_get_lc
  

  function base_get_ctxt(idxmap) result(val)
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer :: val
    
    val = idxmap%ictxt

  end function base_get_ctxt
  

  function base_get_mpic(idxmap) result(val)
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer :: val
    
    val = idxmap%mpic

  end function base_get_mpic
  

  subroutine base_set_state(idxmap,val)
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer, intent(in)  :: val
    
    idxmap%state = val
  end subroutine base_set_state

  subroutine base_set_ctxt(idxmap,val)
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer, intent(in)  :: val
    
    idxmap%ictxt = val
  end subroutine base_set_ctxt
  
  subroutine base_set_gr(idxmap,val)
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer, intent(in)  :: val
    
    idxmap%global_rows = val
  end subroutine base_set_gr

  subroutine base_set_gc(idxmap,val)
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer, intent(in)  :: val
    
    idxmap%global_cols = val
  end subroutine base_set_gc

  subroutine base_set_lr(idxmap,val)
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer, intent(in)  :: val
    
    idxmap%local_rows = val
  end subroutine base_set_lr

  subroutine base_set_lc(idxmap,val)
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer, intent(in)  :: val
    
    idxmap%local_cols = val
  end subroutine base_set_lc

  subroutine base_set_mpic(idxmap,val)
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer, intent(in)  :: val
    
    idxmap%mpic = val
  end subroutine base_set_mpic

  
  function base_row_extendable(idxmap) result(val)
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    logical :: val
    val = .false.
  end function base_row_extendable

  function base_is_repl(idxmap) result(val)
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    logical :: val
    val = .false.
  end function base_is_repl
    
  function base_is_null(idxmap) result(val)
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    logical :: val
    val =  (idxmap%state == psb_desc_null_)
  end function base_is_null
  
  
  function base_is_bld(idxmap) result(val)
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    logical :: val
    val = (idxmap%state == psb_desc_bld_).or.&
         & (idxmap%state == psb_desc_ovl_bld_)
  end function base_is_bld
    
  function base_is_upd(idxmap) result(val)
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    logical :: val
    val = (idxmap%state == psb_desc_upd_)
  end function base_is_upd
    
  function base_is_asb(idxmap) result(val)
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    logical :: val
    val = (idxmap%state == psb_desc_asb_).or.&
         & (idxmap%state == psb_desc_ovl_asb_)
  end function base_is_asb
    
  function base_is_valid(idxmap) result(val)
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    logical :: val
    val = idxmap%is_bld().or.idxmap%is_upd().or.idxmap%is_asb()
  end function base_is_valid

    
  function base_is_ovl(idxmap) result(val)
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    logical :: val
    val = (idxmap%state == psb_desc_ovl_bld_).or.&
         & (idxmap%state == psb_desc_ovl_asb_)
  end function base_is_ovl
  
  function base_sizeof(idxmap) result(val)
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer(psb_long_int_k_) :: val
    
    val = 8 * psb_sizeof_int
  end function base_sizeof


  ! !!!!!!!!!!!!!!!!
  !
  ! !!!!!!!!!!!!!!!!
  subroutine base_l2gs1(idx,idxmap,info,mask,owned)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer, intent(inout) :: idx
    integer, intent(out)   :: info 
    logical, intent(in), optional :: mask
    logical, intent(in), optional :: owned
    Integer :: err_act
    character(len=20)  :: name='base_l2g'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end subroutine base_l2gs1

  subroutine base_l2gs2(idxin,idxout,idxmap,info,mask,owned)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer, intent(in)    :: idxin
    integer, intent(out)   :: idxout
    integer, intent(out)   :: info 
    logical, intent(in), optional :: mask
    logical, intent(in), optional :: owned

    Integer :: err_act
    character(len=20)  :: name='base_l2g'
    logical, parameter :: debug=.false.
    
    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())
    
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return
    
  end subroutine base_l2gs2


  subroutine base_l2gv1(idx,idxmap,info,mask,owned)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer, intent(inout) :: idx(:)
    integer, intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    logical, intent(in), optional :: owned
    Integer :: err_act
    character(len=20)  :: name='base_l2g'
    logical, parameter :: debug=.false.
    
    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())
    
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return
  end subroutine base_l2gv1

  subroutine base_l2gv2(idxin,idxout,idxmap,info,mask,owned)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer, intent(in)    :: idxin(:)
    integer, intent(out)   :: idxout(:)
    integer, intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    logical, intent(in), optional :: owned
    Integer :: err_act
    character(len=20)  :: name='base_l2g'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end subroutine base_l2gv2


  subroutine base_g2ls1(idx,idxmap,info,mask,owned)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer, intent(inout) :: idx
    integer, intent(out)   :: info 
    logical, intent(in), optional :: mask
    logical, intent(in), optional :: owned
    Integer :: err_act
    character(len=20)  :: name='base_g2l'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end subroutine base_g2ls1

  subroutine base_g2ls2(idxin,idxout,idxmap,info,mask,owned)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer, intent(in)    :: idxin
    integer, intent(out)   :: idxout
    integer, intent(out)   :: info 
    logical, intent(in), optional :: mask
    logical, intent(in), optional :: owned

    Integer :: err_act
    character(len=20)  :: name='base_g2l'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return
    
  end subroutine base_g2ls2


  subroutine base_g2lv1(idx,idxmap,info,mask,owned)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer, intent(inout) :: idx(:)
    integer, intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    logical, intent(in), optional :: owned
    Integer :: err_act
    character(len=20)  :: name='base_g2l'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end subroutine base_g2lv1

  subroutine base_g2lv2(idxin,idxout,idxmap,info,mask,owned)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    integer, intent(in)    :: idxin(:)
    integer, intent(out)   :: idxout(:)
    integer, intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    logical, intent(in), optional :: owned

    Integer :: err_act
    character(len=20)  :: name='base_g2l'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return


  end subroutine base_g2lv2



  subroutine base_g2ls1_ins(idx,idxmap,info,mask)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer, intent(inout) :: idx
    integer, intent(out)   :: info 
    logical, intent(in), optional :: mask
    Integer :: err_act
    character(len=20)  :: name='base_g2l_ins'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end subroutine base_g2ls1_ins

  subroutine base_g2ls2_ins(idxin,idxout,idxmap,info,mask)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer, intent(in)    :: idxin
    integer, intent(out)   :: idxout
    integer, intent(out)   :: info 
    logical, intent(in), optional :: mask
    
    Integer :: err_act
    character(len=20)  :: name='base_g2l_ins'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return
    
  end subroutine base_g2ls2_ins


  subroutine base_g2lv1_ins(idx,idxmap,info,mask)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer, intent(inout) :: idx(:)
    integer, intent(out)   :: info 
    logical, intent(in), optional :: mask(:)

    Integer :: err_act
    character(len=20)  :: name='base_g2l_ins'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end subroutine base_g2lv1_ins

  subroutine base_g2lv2_ins(idxin,idxout,idxmap,info,mask)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer, intent(in)    :: idxin(:)
    integer, intent(out)   :: idxout(:)
    integer, intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    Integer :: err_act
    character(len=20)  :: name='base_g2l_ins'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end subroutine base_g2lv2_ins


  subroutine base_asb(idxmap,info)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer, intent(out) :: info
    
    Integer :: err_act
    character(len=20)  :: name='base_asb'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return
    
  end subroutine base_asb

  subroutine base_free(idxmap)
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    
    Integer :: err_act
    character(len=20)  :: name='base_free'
    logical, parameter :: debug=.false.

    ! almost nothing to be done here
    idxmap%state          = -1 
    idxmap%ictxt          = -1
    idxmap%mpic           = -1
    idxmap%global_rows    = -1
    idxmap%global_cols    = -1
    idxmap%local_rows     = -1
    idxmap%local_cols     = -1

    return

  end subroutine base_free

  subroutine base_set_null(idxmap)
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap

    idxmap%state          = psb_desc_null_
    idxmap%ictxt          = -1
    idxmap%mpic           = -1
    idxmap%global_rows    = -1
    idxmap%global_cols    = -1
    idxmap%local_rows     = -1
    idxmap%local_cols     = -1

  end subroutine base_set_null


  function base_get_fmt(idxmap) result(res)
    implicit none 
    class(psb_indx_map), intent(in) :: idxmap
    character(len=5) :: res
    res = 'NULL'
  end function base_get_fmt


  subroutine base_init_vl(idxmap,ictxt,vl,info)
    use psb_penv_mod
    use psb_error_mod
    implicit none 
    class(psb_indx_map), intent(inout) :: idxmap
    integer, intent(in)  :: ictxt, vl(:)
    integer, intent(out) :: info
    Integer :: err_act
    character(len=20)  :: name='base_init_vl'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(psb_err_missing_override_method_,&
         & name,a_err=idxmap%get_fmt())

    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return
  end subroutine base_init_vl
    


end module psb_indx_map_mod
