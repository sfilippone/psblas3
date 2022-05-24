module psb_cutil_cbind_mod

   use iso_c_binding
   use psb_util_mod
   use psb_base_mod
   use psb_objhandle_mod
   use psb_base_string_cbind_mod

contains

    function psb_c_cmm_mat_write(ah,matrixtitle,filename) bind(c) result(res)
        use psb_base_mod
        use psb_util_mod
        use psb_base_string_cbind_mod
        implicit none
        integer(psb_c_ipk_) :: res

        type(psb_c_cspmat)   :: ah
        character(c_char)      :: matrixtitle(*)
        character(c_char)      :: filename(*)

        type(psb_cspmat_type), pointer :: ap
        character(len=1024)     :: mtitle
        character(len=1024)     :: fname
        integer(psb_c_ipk_)     :: info


        res = -1
        if (c_associated(ah%item)) then
          call c_f_pointer(ah%item,ap)
        else
          return
        end if

        call stringc2f(matrixtitle,mtitle)
        call stringc2f(filename,fname)

        call mm_mat_write(ap,mtitle,info,filename=fname)

        res = info

    end function psb_c_cmm_mat_write

    function psb_c_cglobal_mat_write(ah,cdh) bind(c) result(res)
        use psb_base_mod
        use psb_util_mod
        use psb_base_string_cbind_mod
        implicit none
        integer(psb_c_ipk_) :: res

        type(psb_c_cspmat)   :: ah 
        type(psb_c_descriptor) :: cdh

	type(psb_cspmat_type), pointer :: ap
	type(psb_desc_type), pointer :: descp
	! Local variables
	type(psb_cspmat_type) :: aglobal
        integer(psb_ipk_) :: info, iam, np
	type(psb_ctxt_type) :: ctxt
	character(len=40)  :: matrixname

	res = -1
    	if (c_associated(cdh%item)) then
      	  call c_f_pointer(cdh%item,descp)
    	else
      	  return
    	end if
        if (c_associated(ah%item)) then
          call c_f_pointer(ah%item,ap)
        else
         return
        end if

	ctxt = descp%get_ctxt()
	call psb_info(ctxt,iam,np)
	call psb_gather(aglobal,ap,descp,info)
        if (iam == psb_root_) then
	   write(matrixname,'("A-np-",I1,".mtx")') np
      	   call mm_mat_write(aglobal,"Global matrix",info,filename=trim(matrixname))
	end if

	call psb_spfree(aglobal,descp,info)
	res = info

    end function psb_c_cglobal_mat_write

  function psb_c_cglobal_vec_write(vh,cdh) bind(c) result(res)
        use psb_base_mod
        use psb_util_mod
        use psb_base_string_cbind_mod
        implicit none
        integer(psb_c_ipk_) :: res

        type(psb_c_cvector)  :: vh 
        type(psb_c_descriptor) :: cdh

	type(psb_c_vect_type), pointer :: vp
	type(psb_desc_type), pointer :: descp
	! Local variables
	complex(psb_spk_), allocatable :: vglobal(:)
        integer(psb_ipk_) :: info, iam, np
	type(psb_ctxt_type) :: ctxt
	character(len=40)  :: vecname

	res = -1
    	if (c_associated(cdh%item)) then
      	  call c_f_pointer(cdh%item,descp)
    	else
      	  return
    	end if
        if (c_associated(vh%item)) then
          call c_f_pointer(vh%item,vp)
        else
         return
        end if

	ctxt = descp%get_ctxt()
	call psb_info(ctxt,iam,np)
	call psb_gather(vglobal,vp,descp,info)
        if (iam == psb_root_) then
	   write(vecname,'("v-np-",I1,".mtx")') np
      	   call mm_array_write(vglobal,"Global vector",info,filename=trim(vecname))
	end if

	deallocate(vglobal,stat=info)
	res = info

    end function psb_c_cglobal_vec_write


end module psb_cutil_cbind_mod
