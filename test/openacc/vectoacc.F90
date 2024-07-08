program vectoacc
  use psb_base_mod
  use psb_oacc_mod
  implicit none

  type(psb_d_vect_oacc) :: v3, v4, v5
  integer(psb_ipk_) :: info, n, i
  real(psb_dpk_) :: alpha, beta, result
  double precision, external :: etime

  real(psb_dpk_) :: dot_host, dot_dev, t_host, t_dev
  double precision :: time_start, time_end
  integer(psb_ipk_), parameter :: ntests=80, ngpu=20

  write(*, *) 'Test of the vector operations with OpenACC'

  write(*, *) 'Enter the size of the vectors'
  read(*, *) n
  alpha = 2.0
  beta = 0.5

  call v3%all(n, info)
  call v4%all(n, info)
  call v5%all(n, info)

  do i = 1, n
      v3%v(i) = real(i, psb_dpk_)
      v4%v(i) = real(n - i, psb_dpk_)
  end do

  call v3%set_dev()
  call v4%set_dev()

  call v3%scal(alpha)
  call v3%sync()

  do i = 1, n
      if (v3%v(i) /= alpha * real(i, psb_dpk_)) then
          write(*, *) 'Scal error : index', i
      end if
  end do
  write(*, *) 'Scal test passed'

  result = v3%dot_v(n, v4)
  call v3%sync()
  call v4%sync()
  if (result /= sum(v3%v * v4%v)) then
      write(*, *) 'Dot_v error,  expected result:', sum(v3%v * v4%v), 'instead of :', result
  end if
  write(*, *) 'Dot_v test passed'

  result = v3%nrm2(n)
  call v3%sync()
  if (result /= sqrt(sum(v3%v ** 2))) then
      write(*, *) 'nrm2 error, expected result:', sqrt(sum(v3%v ** 2)), 'instead of :', result
  end if
  write(*, *) 'nrm2 test passed'

  call v3%set_host()
  call v4%set_host()

  time_start = etime()
  do i = 1, ntests
      dot_host = sum(v3%v * v4%v)
  end do
  time_end = etime()
  t_host = (time_end - time_start) / real(ntests)
  write(*, *) 'Performance host: ', t_host, ' sec'

  call v3%set_dev()
  call v4%set_dev()
  time_start = etime()
  do i = 1, ntests
      dot_dev = v3%dot_v(n, v4)
  end do
  !$acc wait
  time_end = etime()
  t_dev = (time_end - time_start) / real(ntests)
  write(*, *) 'Performance device: ', t_dev, ' sec'

  call v3%free(info)
  call v4%free(info)
  call v5%free(info)

end program vectoacc
