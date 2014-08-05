subroutine adj_read (a,filename,desc_a,info)

    type(psb_dspmat_type), intent (out) :: a
    character(len=20) :: filename
    type(psb_desc_type):: desc_a
    integer (psb_ipk_) :: info
    implicit none

    integer :: i,j,nnzero,nbCol
    integer :: unitFile, iError,line
    integer(psb_ipk_), allocatable :: fileContent(:,:),ia(:),ja(:)
    real (psb_dpk_), allocatable :: val(:)
    integer(psb_ipk_), allocatable :: size_mat(:)
    
    nbCol = 2
    allocate (size_mat(nbCol))
    unitFile = 1
    open(UNIT=unitFile, FILE=filename, FORM="FORMATTED", STATUS="OLD",
ACTION="READ")

    nnzero = size_mat(2)
    allocate (fileContent(nnzeros,nbCol))

    do line = 1,nnzero
        read(unitFile, *) fileContent(line,1:nbCol)
    end do saveNodes
    close(UNIT=unitFile)

    allocate(ia(nnzero),ja(nnzero),val(nnzero))
    do i=1,nnzero
        ia(i)=fileContent(i,1)
        ja(i)=filecontent(i,2)
        val(i)=1.0
    end do

    call psb_spins(nnzero, ia, ja, val, a, desc_a, info)    

end subroutine adj_read
