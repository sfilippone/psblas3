program datavect
    use psb_base_mod
    use psb_oacc_mod
    implicit none

    type(psb_d_vect_oacc) :: v3, v4, v5
    integer(psb_ipk_) :: info, n, i, old_percentage, percentage
    real(psb_dpk_) :: alpha, dot_dev, dot_host, t_alloc_host, t_alloc_dev, t_calc_host, t_calc_dev
    double precision, external :: etime
    double precision :: time_start, time_end
    integer, parameter :: min_size = 1000, max_size = 100000000, step_size = 1000000
    integer, parameter :: ntests = 80, ngpu = 20
    integer :: size
    character(len=20) :: filename
  
    open(unit=10, file='performance_data.csv', status='unknown')
    write(10, '(A, A, A, A, A)') 'Size,Alloc_Host,Alloc_Dev,Calc_Host,Calc_Dev'
  
    write(*, *) 'Test of the vector operations with OpenACC'
  
    alpha = 2.0
    old_percentage = 0

    do size = min_size, max_size, step_size
        n = size
        percentage = int(real(size - min_size) / real(max_size - min_size) * 100.0)
        if (percentage /= old_percentage) then
            write(*, '(A,I3,A)', advance='no') 'Progress: ', percentage, '%'
            write(*,'(A)', advance='no') char(13) 
            old_percentage = percentage
        end if
    
        time_start = etime()
        call v3%all(n, info)
        call v4%all(n, info)
        call v5%all(n, info)
        time_end = etime()
        t_alloc_host = (time_end - time_start)
    
        do i = 1, n
            v3%v(i) = real(i, psb_dpk_)
            v4%v(i) = real(n - i, psb_dpk_)
        end do

        call v3%scal(alpha)
    
        call v3%set_host()
        call v4%set_host()

        time_start = etime()
        do i = 1, ntests
            dot_host = sum(v3%v * v4%v)
        end do
        time_end = etime()
        t_calc_host = (time_end - time_start) / real(ntests)
        
        time_start = etime()
        call v3%set_dev()
        call v4%set_dev()
        call v3%sync_space()
        call v4%sync_space()
        time_end = etime()
        t_alloc_dev = (time_end - time_start)

        time_start = etime()
        do i = 1, ntests
            dot_dev = v3%dot_v(n, v4)
        end do
        !$acc wait
        time_end = etime()
        t_calc_dev = (time_end - time_start) / real(ntests)
    
        write(10, '(I10, 1X, ES12.5, 1X, ES12.5, 1X, ES12.5, 1X, ES12.5)') size, t_alloc_host, t_alloc_dev, t_calc_host, t_calc_dev
    
        call v3%free(info)
        call v4%free(info)
        call v5%free(info)
    end do
    
    close(10)
    write(*, *) 'Performance data written to performance_data.csv'
  

end program datavect
