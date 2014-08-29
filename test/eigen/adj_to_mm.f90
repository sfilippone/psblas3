program adj_to_mm 
    use psb_base_mod

    ! reads a file containing the coordinate description of an adjacence
    ! sparse matrix (values are always 1) and stores it into 'a' 
    !
    ! on entry :
    ! file_adj : file containing the description of the adjacency matrix
    ! file_mm  : destination MM format file containing the same matrix
    ! begin    : integer specifying the first index of rows and columns in the
    !            file_adj (usually 1 or 0)

    character(len=40) :: file_adj,file_mm
    integer(psb_ipk_) :: i,j,k,nnzero,nrows,begin
    type(psb_d_coo_sparse_mat) :: acoo

    read(psb_inp_unit,*)file_adj
    read(psb_inp_unit,*)file_mm
    read(psb_inp_unit,*)begin

    open(15, FILE=file_adj, STATUS="OLD", ACTION="READ")
    open(14, FILE=file_mm, ACTION="WRITE")
    read(15, *) nrows,nnzero
    write(14,'(i20,i20,i20)') nrows,nrows,nnzero

    do k = 1,nnzero
        read(15, *) i,j
        j=j+(1-begin)
        i=i+(1-begin)
        write(14,'(i20,i20,F20.3)')i,j,1.0
    end do
    close(UNIT=15)
    close(UNIT=14)
end program adj_to_mm

