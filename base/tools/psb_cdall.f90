subroutine psb_cdall(ictxt, desc, info,mg,ng,parts,vg,vl,flag,nl,repl,globalcheck,lidx,usehash)
  use psb_desc_mod
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod
  use psb_cd_tools_mod, psb_protect_name => psb_cdall
  use psi_mod
  implicit None
  procedure(psb_parts)              :: parts
  integer(psb_ipk_), intent(in)     :: mg,ng,ictxt, vg(:), vl(:),nl,lidx(:)
  integer(psb_ipk_), intent(in)     :: flag
  logical, intent(in)               :: repl, globalcheck,usehash
  integer(psb_ipk_), intent(out)    :: info
  type(psb_desc_type), intent(out)  :: desc

  optional :: mg,ng,parts,vg,vl,flag,nl,repl, globalcheck,lidx, usehash

  interface 
    subroutine psb_cdals(m, n, parts, ictxt, desc, info)
      use psb_desc_mod
      procedure(psb_parts)              :: parts
      integer(psb_ipk_), intent(in)     :: m,n,ictxt
      Type(psb_desc_type), intent(out)  :: desc
      integer(psb_ipk_), intent(out)    :: info
    end subroutine psb_cdals
    subroutine psb_cdalv(v, ictxt, desc, info, flag)
      use psb_desc_mod
      integer(psb_ipk_), intent(in)           :: ictxt, v(:)
      integer(psb_ipk_), intent(in), optional :: flag
      integer(psb_ipk_), intent(out)          :: info
      Type(psb_desc_type), intent(out)        :: desc
    end subroutine psb_cdalv
    subroutine psb_cd_inloc(v, ictxt, desc, info, globalcheck,idx, usehash)
      use psb_desc_mod
      implicit None
      integer(psb_ipk_), intent(in)           :: ictxt, v(:)
      integer(psb_ipk_), intent(out)          :: info
      type(psb_desc_type), intent(out)        :: desc
      logical, intent(in), optional           :: globalcheck, usehash
      integer(psb_ipk_), intent(in), optional :: idx(:)
    end subroutine psb_cd_inloc
    subroutine psb_cdrep(m, ictxt, desc,info)
      use psb_desc_mod
      integer(psb_ipk_), intent(in)     :: m,ictxt
      Type(psb_desc_type), intent(out)  :: desc
      integer(psb_ipk_), intent(out)    :: info
    end subroutine psb_cdrep
  end interface
  character(len=20)  :: name
  integer(psb_ipk_)  :: err_act, n_, flag_, i, me, np, nlp, nnv, lr
  logical            :: usehash_
  integer(psb_ipk_), allocatable :: itmpsz(:) 
  integer(psb_mpik_) :: iictxt
 
  

  if (psb_get_errstatus() /= 0) return 
  info=psb_success_
  name = 'psb_cdall'
  call psb_erractionsave(err_act)

  call psb_info(ictxt, me, np)
  iictxt  = ictxt
  if (count((/ present(vg),present(vl),&
       &  present(parts),present(nl), present(repl) /)) /= 1) then 
    info=psb_err_no_optional_arg_
    call psb_errpush(info,name,a_err=" vg, vl, parts, nl, repl")
    goto 9999 
  endif

  desc%base_desc => null() 
  if (allocated(desc%indxmap)) then 
    write(0,*) 'Allocated on an intent(OUT) var?'
  end if

  if (present(parts)) then 

    if (.not.present(mg)) then 
      info=psb_err_no_optional_arg_
      call psb_errpush(info,name)
      goto 9999 
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
      goto 9999 
    end if
    if (.not.repl) then 
      info=psb_err_no_optional_arg_
      call psb_errpush(info,name)
      goto 9999 
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

    call psb_cd_inloc(vl(1:nnv),ictxt,desc,info, globalcheck=globalcheck,idx=lidx)

  else if (present(nl)) then 
    
    if (present(usehash)) then
      usehash_ = usehash
    else
      usehash_ = .false.
    end if
    if (usehash_) then
      write(0,*) 'Fix usehash_ implementationt '
    end if

    if (np == 1) then 
      allocate(psb_repl_map      :: desc%indxmap, stat=info)
    else
      allocate(psb_gen_block_map :: desc%indxmap, stat=info)
    end if
    if (info == psb_success_) then 
      select type(aa => desc%indxmap) 
      type is (psb_repl_map) 
        call aa%repl_map_init(iictxt,nl,info)
      type is (psb_gen_block_map) 
        call aa%gen_block_map_init(iictxt,nl,info)
      class default 
        ! This cannot happen 
        info = psb_err_internal_error_
        goto 9999
      end select
    end if
    
    call psb_realloc(1,itmpsz, info)
    if (info /= 0) then 
      write(0,*) 'Error reallocating itmspz'
      goto 9999
    end if
    itmpsz(:) = -1
    call psi_bld_tmpovrl(itmpsz,desc,info)

  endif

  if (info /= psb_success_) goto 9999

  ! Finish off 
  lr = desc%indxmap%get_lr()
  call psb_realloc(max(1,lr/2),desc%halo_index, info)
  if (info == psb_success_) call psb_realloc(1,desc%ext_index, info)
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='psb_realloc')
    Goto 9999
  end if
  desc%halo_index(:)           = -1
  desc%ext_index(:)            = -1
  call psb_cd_set_bld(desc,info)
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psb_cdall
