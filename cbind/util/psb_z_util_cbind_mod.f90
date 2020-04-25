module psb_zutil_cbind_mod

   use iso_c_binding
   use psb_util_mod
   use psb_base_mod
   use psb_objhandle_mod
   use psb_base_string_cbind_mod

contains

    function psb_c_zmm_mat_write(ah,matrixtitle,filename) bind(c) result(res)
        use psb_base_mod
        use psb_util_mod
        use psb_base_string_cbind_mod
        implicit none
        integer(psb_c_ipk_) :: res

        type(psb_c_zspmat)   :: ah
        character(c_char)      :: matrixtitle(*)
        character(c_char)      :: filename(*)

        type(psb_zspmat_type), pointer :: ap
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

    end function psb_c_zmm_mat_write

end module psb_zutil_cbind_mod
