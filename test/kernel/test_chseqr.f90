program test_chseqr
        
complex , dimension (1:3,1:3) :: H,Z,Vr
complex , dimension (1:6,1:6) :: Rwork
complex, dimension (1:3) :: eing,work
integer :: info, N
N=3
H = reshape((/ (1.0,0),(0,0), (2,0), (2,0) , (3,0) , (-4,0) , (0,0) , (0,0) , (2,0) /), shape(H))
do i=1,N
        do j=1,N
                write(*,'("H : "i5,i5,2x,g20.4,g20.4)')i,j,real(H(i,j)),aimag(H(i,j))
        end do
end do
call CHSEQR('E','N',N,1,N,H,N,eing,Z,N,work,N,info)
!call CGEEV('N','N',N,H,N,eing,Z,N,Vr,N,work,3*N,Rwork,info)
do i=1,N
       write(*,'("valeur propre de H : "i5, i5, g20.4,g20.4)')info, i,real(eing(i)),aimag(eing(i))
enddo

end program test_chseqr
