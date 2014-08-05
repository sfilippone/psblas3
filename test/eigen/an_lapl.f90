program analyse
        use psb_base_mod
        use psb_util_mod
        implicit none

        character(len=30) :: lapl_file,max_file,res_file,analyse_file,mat_file
        integer (psb_ipk_) :: i,nb_mat, n
        real(psb_dpk_) ::tab(18),eig1,eig2,eig0,eigmin,eigmax,delta,delta2,&
                &deltamax,timemin,timemax,t1,t2,time0,t0

        read(psb_inp_unit,*) lapl_file
        read(psb_inp_unit,*) max_file
        read(psb_inp_unit,*) res_file
        read(psb_inp_unit,*) analyse_file
        read(psb_inp_unit,*) nb_mat

        open(15,FILE=lapl_file,action="read")
        open(13,FILE=max_file,action="read")
        open(14,file=res_file,action="read")

        do i=1,21
                read(15,*) 
        end do
        do i=1,22
                read(13,*)
        end do
        
        open (16, file=analyse_file,position="append",action="write")
        n=0
        delta=0
        delta2=0
        deltamax=0
        timemax=0
        timemin=0
        time0=0
        do i=1,nb_mat
                read(15,*)tab(1:18)
                read(14,*)mat_file,eig1,eig2,t1     
                read(14,*)eigmin,eig0,t2          
                
                if(abs(tab(6)-eigmin)>0.001)then
                         n=n+1
                         write(16,'("smallest of : ",F20.2,F20.2)')tab(2),tab(3)
                         write(16,'(F20.7,F20.7)')tab(6),eigmin
                end if      
                if(abs(tab(16)-eig2)>0.01)then
                         n=n+1
                        write(16,'("lambda N-1 of : ",F20.2,F20.2)')tab(2),tab(3)
                        write(16,'(F20.7,F20.7)')tab(16),eig2
                end if
                if(abs(tab(18)-eig1)>0.01)then
                         n=n+1
                         write(16,'("biggest of : ", F20.2,F20.2)')tab(2),tab(3)
                         write(16,'(F20.7,F20.7)')tab(18),eig1
                end if
                timemax=timemax+t1
                timemin=timemin+t2
                delta2=delta2+abs(1-eigmin/tab(6))
                delta=delta+abs(1-eig2/tab(16))+abs(1-eig1/tab(18))

                read(13,*)tab(1:7)
                read(14,*)mat_file,eigmax,t0
                if(abs(tab(4)-eigmax)>0.01)then
                         write(16,'(F20.2,F20.2)')tab(2),tab(3)
                         write(16,'(F20.6,F20.6)')tab(4),eigmax
                endif
                time0=time0+t0
                deltamax=deltamax+abs(1-eigmax/tab(4))
        end do
        write(16,'("errors on ",i20," eigenvalues")')n
        write(16,'("gap average for lapl max",g20.5)')delta/(2*nb_mat)
        write(16,'("gap average for lapl min",g20.5)')delta2/(nb_mat)
        write(16,'("gap average for max",g20.5)')deltamax/(nb_mat)
        write(16,'("time average lapl max",g20.5)')timemax/(nb_mat)
        write(16,'("time average lapl min",g20.5)')timemin/(nb_mat)
        write(16,'("time average A max",g20.5)')time0/(nb_mat)
        close(15)
        close(14)
        close(13)
        close(16)

end program analyse
