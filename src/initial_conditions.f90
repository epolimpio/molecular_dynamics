module initial_conditions

	implicit none
	private init_random_seed, max_boltz
	public initial_pos, initial_mom


contains

	function initial_mom(N, temp) result(mom)

		real(8), intent(in) :: temp
		integer(8), intent(in) :: N
		real(8) :: rnd



	end function initial_mom


	function max_boltz(N, temp) result(velocity)

		real(8), intent(in) :: temp
		real(8), parameter :: PI = 4*atan(1d0)
		real(8) :: velocity(3, N)
		integer :: i, j
		real(8) :: v_sum, v_rnd, var_v

		! We use Gaussian random to generate the particle velocity 
		! in each direction for N-1 particles. The last particle 
		! velocity is such that the sum is zero
		do i = 1, 3
			v_sum = 0
			do j = 1, N-1
				v_rnd = box_muller
				velocity(i, j) = v_rnd
				v_sum = v_sum + v_rnd
			end do
		end do

	end function max_boltz

	! Implementation of polar Box-Muller algorithm
	! to generate Gaussian Random numbers

	real(8) function box_muller() result(gauss_number)

		real(8) :: x1, x2, w, y1, y2

		call init_random_seed

		do
			x1 = 2d0*random_number - 1
			x2 = 2d0*random_number - 1
			w = x1 * x1 + x2 * x2
		until w <= 1d0

		gauss_number = sqrt( (-2d0 * log(w)) / w )

	end function box_muller


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