subroutine psb_cdall(ictxt, desc, info,mg,ng,parts,vg,vl,flag,nl,repl, globalcheck)
  use psb_descriptor_type
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod
  use psb_base_tools_mod, psb_protect_name => psb_cdall
  use psi_mod
  implicit None
  procedure(psb_parts)               :: parts
  integer(psb_ipk_), intent(in)               :: mg,ng,ictxt, vg(:), vl(:),nl
  integer(psb_ipk_), intent(in)               :: flag
  logical, intent(in)               :: repl, globalcheck
  integer(psb_ipk_), intent(out)              :: info
  type(psb_desc_type), intent(out)  :: desc

  optional :: mg,ng,parts,vg,vl,flag,nl,repl, globalcheck

  interface 
    subroutine psb_cdals(m, n, parts, ictxt, desc, info)
      use psb_descriptor_type
      procedure(psb_parts)               :: parts
      integer(psb_ipk_), intent(in)                 :: m,n,ictxt
      Type(psb_desc_type), intent(out)    :: desc
      integer(psb_ipk_), intent(out)                :: info
    end subroutine psb_cdals
    subroutine psb_cdalv(v, ictxt, desc, info, flag)
      use psb_descriptor_type
      integer(psb_ipk_), intent(in)               :: ictxt, v(:)
      integer(psb_ipk_), intent(in), optional     :: flag
      integer(psb_ipk_), intent(out)              :: info
      Type(psb_desc_type), intent(out)  :: desc
    end subroutine psb_cdalv
    subroutine psb_cd_inloc(v, ictxt, desc, info, globalcheck)
      use psb_descriptor_type
      implicit None
      integer(psb_ipk_), intent(in)               :: ictxt, v(:)
      integer(psb_ipk_), intent(out)              :: info
      type(psb_desc_type), intent(out)  :: desc
      logical, intent(in), optional     :: globalcheck
    end subroutine psb_cd_inloc
    subroutine psb_cdrep(m, ictxt, desc,info)
      use psb_descriptor_type
      integer(psb_ipk_), intent(in)               :: m,ictxt
      Type(psb_desc_type), intent(out)  :: desc
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psb_cdrep
  end interface
  character(len=20)   :: name
  integer(psb_ipk_) :: err_act, n_, flag_, i, me, np, nlp, nnv, lr
  integer(psb_ipk_), allocatable :: itmpsz(:) 



  if (psb_get_errstatus() /= 0) return 
  info=psb_success_
  name = 'psb_cdall'
  call psb_erractionsave(err_act)

  call psb_info(ictxt, me, np)

  if (count((/ present(vg),present(vl),&
       &  present(parts),present(nl), present(repl) /)) /= 1) then 
    info=psb_err_no_optional_arg_
    call psb_errpush(info,name,a_err=" vg, vl, parts, nl, repl")
    goto 999 
  endif

  desc%base_desc => null() 
  if (allocated(desc%indxmap)) then 
    write(0,*) 'Allocated on an intent(OUT) var?'
  end if

  if (present(parts)) then 

    if (.not.present(mg)) then 
      info=psb_err_no_optional_arg_
      call psb_errpush(info,name)
      goto 999 
    end if
    if (present(ng)) then 
      n_ = ng
    else
      n_ = mg 
    endif
    call  psb_cdals(mg, n_, parts, ictxt, desc, info)

  else if (present(repl)) then 

    if (.not.present(mg)) then 
      info=psb_err_no_optional_arg_
      call psb_errpush(info,name)
      goto 999 
    end if
    if (.not.repl) then 
      info=psb_err_no_optional_arg_
      call psb_errpush(info,name)
      goto 999 
    end if

    call  psb_cdrep(mg, ictxt, desc, info)


  else if (present(vg)) then 

    if (present(flag)) then 
      flag_=flag
    else
      flag_=0
    endif
    if (present(mg)) then 
      nnv = min(mg,size(vg))
    else
      nnv = size(vg)
    end if

    call psb_cdalv(vg(1:nnv), ictxt, desc, info, flag=flag_)

  else if (present(vl)) then 

    if (present(nl)) then 
      nnv = min(nl,size(vl))
    else
      nnv = size(vl)
    end if

    call psb_cd_inloc(vl(1:nnv),ictxt,desc,info, globalcheck=globalcheck)

  else if (present(nl)) then 
    

    if (np == 1) then 
      allocate(psb_repl_map      :: desc%indxmap, stat=info)
    else
      allocate(psb_gen_block_map :: desc%indxmap, stat=info)
    end if
    if (info == psb_success_) then 
      select type(aa => desc%indxmap) 
      type is (psb_repl_map) 
        call aa%repl_map_init(ictxt,nl,info)
      type is (psb_gen_block_map) 
        call aa%gen_block_map_init(ictxt,nl,info)
      class default 
        ! This cannot happen 
        info = psb_err_internal_error_
        goto 999
      end select
    end if
    
    call psb_realloc(1,itmpsz, info)
    if (info /= 0) then 
      write(0,*) 'Error reallocating itmspz'
      goto 999
    end if
    itmpsz(:) = -1
    call psi_bld_tmpovrl(itmpsz,desc,info)

  endif

  if (info /= psb_success_) goto 999

  ! Finish off 
  lr = desc%indxmap%get_lr()
  call psb_realloc(max(1,lr/2),desc%halo_index, info)
  if (info == psb_success_) call psb_realloc(1,desc%ext_index, info)
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='psb_realloc')
    Goto 999
  end if
  desc%halo_index(:)           = -1
  desc%ext_index(:)            = -1
  call psb_cd_set_bld(desc,info)
  if (info /= psb_success_) goto 999

  call psb_erractionrestore(err_act)
  return

999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return


end subroutine psb_cdall
