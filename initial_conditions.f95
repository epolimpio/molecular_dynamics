module initial_conditions

  implicit none
  private init_random_seed
  public max_boltz, fcc_lattice

contains
  
  subroutine fcc_lattice(n, r, position)
    
    real(8), intent(in) :: r
    integer, intent(in) :: n
    real(8), intent(out) :: position(3, n)

    integer :: num_cubes, i, j
    real(8) :: cube_corner_pos(3)

    num_cubes = n/4

    position(:,1) = (/ 0, 0, 0 /)
    cube_corner_pos = position(:, 1)

    do i = 2,   

  

  end subroutine lattice

  subroutine max_boltz(n, temp, momenta)
  ! We use Gaussian random to generate the particle momenta 
  ! in each direction for N-1 particles. The last particle 
  ! momenta is such that the sum is zero

    real(8), intent(in) :: temp
    integer, intent(in) :: n
    real(8), intent(out) :: momenta(3, n)
    integer :: i, j
    real(8) :: p_sum, p_rnd, var_p
	  
    var_p = temp

    do i = 1, 3
      p_sum = 0
      do j = 1, n
        call box_muller(p_rnd)
        momenta(i, j) = var_p*p_rnd
      end do
      momenta(i,:) = momenta(i,:) - sum(momenta(i,:))/n
      print *, sum(momenta(i,:))/n
     end do
 	  
  end subroutine max_boltz

  subroutine box_muller(gauss_number)
	! Implementation of polar Box-Muller algorithm
	! to generate Gaussian Random numbers

    real(8) :: x(2)
    real(8), parameter :: TWOPI = 8*atan(1d0)
    real(8), intent(out) :: gauss_number

    call init_random_seed
    call random_number(x)

    gauss_number = sqrt(-2*log(x(1))) * cos(TWOPI * x(2))

  end subroutine box_muller

	! Function to initialize the seed of the random
	! number generator

  subroutine init_random_seed()

    integer, allocatable :: seed(:)
    integer :: i, n, un, istat, dt(8), pid, t(2), s
    integer(8) :: count, tms

    call random_seed(size = n)
    allocate(seed(n))
    open(newunit=un, file="/dev/urandom", access="stream",&
     form="unformatted", action="read", status="old", &
     iostat=istat)
    if (istat == 0) then
      read(un) seed
      close(un)
    else
      call system_clock(count)
      if (count /= 0) then
      t = transfer(count, t)
    else
      call date_and_time(values=dt)
      tms = (dt(1) - 1970)*365_8 * 24 * 60 * 60 * 1000 &
       + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
       + dt(3) * 24 * 60 * 60 * 60 * 1000 &
       + dt(5) * 60 * 60 * 1000 &
       + dt(6) * 60 * 1000 + dt(7) * 1000 &
       + dt(8)
      t = transfer(tms, t)
    end if
      s = ieor(t(1), t(2))
      pid = getpid() + 1099279 ! Add a prime
      s = ieor(s, pid)
      if (n >= 3) then
        seed(1) = t(1) + 36269
        seed(2) = t(2) + 72551
        seed(3) = pid
        if (n > 3) then
          seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
        end if
      else
        seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
      end if
    end if
    call random_seed(put=seed)

  end subroutine init_random_seed

end module
