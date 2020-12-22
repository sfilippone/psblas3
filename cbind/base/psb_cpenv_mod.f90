module psb_cpenv_mod
  use iso_c_binding
  use psb_objhandle_mod

  integer, private :: psb_c_index_base=0

contains

  function psb_c_get_index_base() bind(c) result(res)
    implicit none

    integer(psb_c_ipk_)  :: res

    res = psb_c_index_base
  end function psb_c_get_index_base

  subroutine psb_c_set_index_base(base) bind(c)
    implicit none
    integer(psb_c_ipk_), value  :: base

    psb_c_index_base = base
  end subroutine psb_c_set_index_base

  function psb_c_get_errstatus() bind(c) result(res)
    use psb_base_mod, only : psb_get_errstatus, psb_ctxt_type
    implicit none

    integer(psb_c_ipk_)  :: res

    res = psb_get_errstatus()
  end function psb_c_get_errstatus

  subroutine psb_c_init(cctxt) bind(c)
    use psb_base_mod, only : psb_init, psb_ctxt_type
    implicit none

    type(psb_c_object_type)      :: cctxt
    type(psb_ctxt_type), pointer :: ctxt
    integer :: info

    if (c_associated(cctxt%item)) then
      call c_f_pointer(cctxt%item,ctxt)
      deallocate(ctxt,stat=info)
      if (info /= 0) return
    end if
    allocate(ctxt,stat=info)
    if (info /= 0) return
    call psb_init(ctxt)
    cctxt%item = c_loc(ctxt)

  end subroutine psb_c_init

  function psb_c2f_ctxt(cctxt)   result(res)
    implicit none
    type(psb_c_object_type), value :: cctxt
    type(psb_ctxt_type),   pointer :: res

    !res%ctxt = cctxt%ctxt
    if  (.not.c_associated(cctxt%item)) then
      write(0,*) 'Null item in c2f_ctxt? '
      flush(0)
    end if
    if (c_associated(cctxt%item)) call c_f_pointer(cctxt%item,res)
  end function psb_c2f_ctxt

  subroutine psb_c_get_i_ctxt(cctxt,ictxt,info) bind(c)
    implicit none
    type(psb_c_object_type), value :: cctxt
    integer(psb_c_ipk_)  :: ictxt
    integer(psb_c_ipk_)  :: info

    ! Local variables
    type(psb_ctxt_type),   pointer :: ctxt

    ctxt => psb_c2f_ctxt(cctxt)

    call ctxt%get_i_ctxt(ictxt,info)

  end subroutine

  function psb_c_cmp_ctxt(cctxt1, cctxt2) bind(c,name="psb_c_cmp_ctxt") result(res)
    implicit none
    type(psb_c_object_type), value :: cctxt1, cctxt2
    logical(c_bool) :: res

    logical :: equal
    type(psb_ctxt_type),   pointer :: ctxt1, ctxt2

    ctxt1 => psb_c2f_ctxt(cctxt1)
    ctxt2 => psb_c2f_ctxt(cctxt2)

    equal = psb_cmp_ctxt(ctxt1, ctxt2)

    res = equal

  end function

  subroutine psb_c_exit_ctxt(cctxt) bind(c)
    use psb_base_mod, only : psb_exit, psb_ctxt_type
    type(psb_c_object_type), value :: cctxt

    type(psb_ctxt_type), pointer :: ctxt
    ctxt => psb_c2f_ctxt(cctxt)
    call psb_exit(ctxt,close=.false.)
    return
  end subroutine psb_c_exit_ctxt

  subroutine psb_c_exit(cctxt) bind(c)
    use psb_base_mod, only : psb_exit, psb_ctxt_type
    type(psb_c_object_type), value :: cctxt

    type(psb_ctxt_type), pointer :: ctxt
    ctxt => psb_c2f_ctxt(cctxt)

    call psb_exit(ctxt)
    return
  end subroutine psb_c_exit

  subroutine psb_c_abort(cctxt) bind(c)
    use psb_base_mod, only : psb_abort, psb_ctxt_type
    type(psb_c_object_type), value :: cctxt

    type(psb_ctxt_type), pointer :: ctxt
    ctxt => psb_c2f_ctxt(cctxt)
    call psb_abort(ctxt)
    return
  end subroutine psb_c_abort


  subroutine psb_c_info(cctxt,iam,np) bind(c)
    use psb_base_mod, only : psb_info, psb_ctxt_type
    type(psb_c_object_type), value :: cctxt
    integer(psb_c_ipk_)        :: iam,np

    type(psb_ctxt_type), pointer :: ctxt
    ctxt => psb_c2f_ctxt(cctxt)
    call psb_info(ctxt,iam,np)
    return
  end subroutine psb_c_info

  subroutine psb_c_barrier(cctxt) bind(c)
    use psb_base_mod, only : psb_barrier, psb_ctxt_type
    type(psb_c_object_type), value :: cctxt

    type(psb_ctxt_type), pointer :: ctxt
    ctxt => psb_c2f_ctxt(cctxt)
    call psb_barrier(ctxt)
  end subroutine psb_c_barrier

  real(c_double) function psb_c_wtime() bind(c)
    use psb_base_mod, only : psb_wtime, psb_ctxt_type

    psb_c_wtime = psb_wtime()
  end function psb_c_wtime

  subroutine psb_c_mbcast(cctxt,n,v,root) bind(c)
    use psb_base_mod, only : psb_bcast, psb_ctxt_type
    implicit none
    type(psb_c_object_type), value :: cctxt
    integer(psb_c_ipk_), value     :: n, root
    integer(psb_c_mpk_)        :: v(*)

    type(psb_ctxt_type), pointer :: ctxt
    ctxt => psb_c2f_ctxt(cctxt)

    if (n < 0) then
      write(0,*) 'Wrong size in BCAST'
      return
    end if
    if (n==0) return

    call psb_bcast(ctxt,v(1:n),root=root)
  end subroutine psb_c_mbcast

  subroutine psb_c_ibcast(cctxt,n,v,root) bind(c)
    use psb_base_mod, only : psb_bcast, psb_ctxt_type
    implicit none
    type(psb_c_object_type), value :: cctxt
    integer(psb_c_ipk_), value     :: n, root
    integer(psb_c_ipk_)        :: v(*)
    type(psb_ctxt_type), pointer :: ctxt

    ctxt => psb_c2f_ctxt(cctxt)

    if (n < 0) then
      write(0,*) 'Wrong size in BCAST'
      return
    end if
    if (n==0) return

    call psb_bcast(ctxt,v(1:n),root=root)
  end subroutine psb_c_ibcast

  subroutine psb_c_lbcast(cctxt,n,v,root) bind(c)
    use psb_base_mod, only : psb_bcast, psb_ctxt_type
    implicit none
    type(psb_c_object_type), value :: cctxt
    integer(psb_c_ipk_), value     :: n, root
    integer(psb_c_lpk_)        :: v(*)
    type(psb_ctxt_type), pointer :: ctxt
    ctxt => psb_c2f_ctxt(cctxt)

    if (n < 0) then
      write(0,*) 'Wrong size in BCAST'
      return
    end if
    if (n==0) return

    call psb_bcast(ctxt,v(1:n),root=root)
  end subroutine psb_c_lbcast

  subroutine psb_c_ebcast(cctxt,n,v,root) bind(c)
    use psb_base_mod, only : psb_bcast, psb_ctxt_type
    implicit none
    type(psb_c_object_type), value :: cctxt
    integer(psb_c_ipk_), value     :: n, root
    integer(psb_c_epk_)        :: v(*)
    type(psb_ctxt_type), pointer :: ctxt
    ctxt => psb_c2f_ctxt(cctxt)

    if (n < 0) then
      write(0,*) 'Wrong size in BCAST'
      return
    end if
    if (n==0) return

    call psb_bcast(ctxt,v(1:n),root=root)
  end subroutine psb_c_ebcast

  subroutine psb_c_sbcast(cctxt,n,v,root) bind(c)
    use psb_base_mod
    implicit none
    type(psb_c_object_type), value :: cctxt
    integer(psb_c_ipk_), value     :: n, root
    real(c_float)     :: v(*)
    type(psb_ctxt_type), pointer :: ctxt
    ctxt => psb_c2f_ctxt(cctxt)

    if (n < 0) then
      write(0,*) 'Wrong size in BCAST'
      return
    end if
    if (n==0) return

    call psb_bcast(ctxt,v(1:n),root=root)
  end subroutine psb_c_sbcast

  subroutine psb_c_dbcast(cctxt,n,v,root) bind(c)
    use psb_base_mod, only : psb_bcast, psb_ctxt_type
    implicit none
    type(psb_c_object_type), value :: cctxt
    integer(psb_c_ipk_), value     :: n, root
    real(c_double)        :: v(*)
    type(psb_ctxt_type), pointer :: ctxt
    ctxt => psb_c2f_ctxt(cctxt)

    if (n < 0) then
      write(0,*) 'Wrong size in BCAST'
      return
    end if
    if (n==0) return

    call psb_bcast(ctxt,v(1:n),root=root)
  end subroutine psb_c_dbcast


  subroutine psb_c_cbcast(cctxt,n,v,root) bind(c)
    use psb_base_mod, only : psb_bcast, psb_ctxt_type
    implicit none
    type(psb_c_object_type), value :: cctxt
    integer(psb_c_ipk_), value     :: n, root
    complex(c_float_complex)        :: v(*)
    type(psb_ctxt_type), pointer :: ctxt
    ctxt => psb_c2f_ctxt(cctxt)

    if (n < 0) then
      write(0,*) 'Wrong size in BCAST'
      return
    end if
    if (n==0) return

    call psb_bcast(ctxt,v(1:n),root=root)
  end subroutine psb_c_cbcast

  subroutine psb_c_zbcast(cctxt,n,v,root) bind(c)
    use psb_base_mod
    implicit none
    type(psb_c_object_type), value :: cctxt
    integer(psb_c_ipk_), value     :: n, root
    complex(c_double_complex)     :: v(*)
    type(psb_ctxt_type), pointer :: ctxt
    ctxt => psb_c2f_ctxt(cctxt)

    if (n < 0) then
      write(0,*) 'Wrong size in BCAST'
      return
    end if
    if (n==0) return

    call psb_bcast(ctxt,v(1:n),root=root)
  end subroutine psb_c_zbcast

  subroutine psb_c_hbcast(cctxt,v,root) bind(c)
    use psb_base_mod, only : psb_bcast, psb_info, psb_ipk_, psb_ctxt_type
    implicit none
    type(psb_c_object_type), value :: cctxt
    integer(psb_c_ipk_), value     :: root
    character(c_char)  :: v(*)
    integer(psb_ipk_)  :: iam, np, n
    type(psb_ctxt_type), pointer :: ctxt
    ctxt => psb_c2f_ctxt(cctxt)

    call psb_info(ctxt,iam,np)

    if (iam==root) then
      n = 1
      do
        if (v(n) == c_null_char) exit
        n = n + 1
      end do
    end if
    call psb_bcast(ctxt,n,root=root)
    call psb_bcast(ctxt,v(1:n),root=root)
  end subroutine psb_c_hbcast

  function psb_c_f2c_errmsg(cmesg,len) bind(c) result(res)
    use psb_base_mod, only : psb_errpop,psb_max_errmsg_len_, psb_ctxt_type
    use psb_base_string_cbind_mod
    implicit none
    character(c_char), intent(inout)  :: cmesg(*)
    integer(psb_c_ipk_), intent(in), value :: len
    integer(psb_c_ipk_) :: res
    character(len=psb_max_errmsg_len_), allocatable :: fmesg(:)
    character(len=psb_max_errmsg_len_) :: tmp
    integer :: i, j, ll, il

    res = 0
    call psb_errpop(fmesg)
    ll = 1
    if (allocated(fmesg)) then
      res = size(fmesg)
      do i=1, size(fmesg)
        tmp = fmesg(i)
        il = len_trim(tmp)
        il = min(il,len-ll)
        !write(0,*) 'loop f2c_errmsg: ', ll,il
        call stringf2c(tmp(1:il),cmesg(ll:ll+il))
        cmesg(ll+il)=c_new_line
        ll = ll+il+1
      end do
      !write(0,*) 'From f2c_errmsg: ', ll,len
    end if
    cmesg(ll) = c_null_char
  end function psb_c_f2c_errmsg

  subroutine psb_c_seterraction_ret() bind(c)
    use psb_base_mod, only : psb_set_erraction, psb_act_ret_, psb_ctxt_type
    call psb_set_erraction(psb_act_ret_)
  end subroutine psb_c_seterraction_ret

  subroutine psb_c_seterraction_print() bind(c)
    use psb_base_mod, only : psb_set_erraction, psb_act_print_, psb_ctxt_type
    call psb_set_erraction(psb_act_print_)
  end subroutine psb_c_seterraction_print

  subroutine psb_c_seterraction_abort() bind(c)
    use psb_base_mod, only : psb_set_erraction, psb_act_abort_, psb_ctxt_type
    call psb_set_erraction(psb_act_abort_)
  end subroutine psb_c_seterraction_abort

end module psb_cpenv_mod
