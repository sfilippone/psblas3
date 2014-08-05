program analyse
        use psb_base_mod
        use psb_util_mod
        implicit none

        character(len=30) :: given_file,res_file,analyse_file,mat_file
        integer (psb_ipk_) :: i,nb_mat, n
        real(psb_dpk_) :: tab(7),eig

        read(psb_inp_unit,*) given_file
        read(psb_inp_unit,*) res_file
        read(psb_inp_unit,*) analyse_file
        read(psb_inp_unit,*) nb_mat
        
        
        open(15,FILE=given_file,action="read")
        open(14,file=res_file,action="read")

        do i=1,22
                read(15,*) 
        end do
        
        open (16, file=analyse_file,position="append",action="write")
        n=0
        do i=1,nb_mat
                read(15,*)tab(1:7)
                read(14,*)mat_file,eig                
                      
                if(tab(4)-eig==0)then
                         write(psb_out_unit,'("on a gagne !")')
                         n=n+1
                else 
                         write(psb_out_unit,'("perdu")')
                         write(16,'(F20.2,F20.2)')tab(2),tab(3)
                         write(16,'(F20.7,F20.7)')tab(4),eig
                end if
        end do
        write(psb_out_unit,'(i20)')n
        close(14)
        close(15)
        close(16)

end program analyse
