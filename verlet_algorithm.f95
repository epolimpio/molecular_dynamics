program ArgonSim
use main
implicit none


subroutine verlet_algorithm(gas,N)
	real (kind=ikind), intent(inout), dimension(N,3,3) :: gas
	real (kind=ikind), dimension(N,3,3) :: gas2
	integer, intent(in) :: dt,m,N
	integer dimension(3) :: r,v,a
	do n = 1,N
		do i = 1,3
			r(i)=gas(n,1,i)
			v(i)=gas(n,2,i)/m
			a(i)=gas(n,3,i)/m
		end do

		do i = 1,3
			r(i) = r(i)+v(i)*dt+.5*a(i)*dt*dt
		end do
		do i = 1,3
			v(i) = v(i)+.5*a(i)*dt !v(t+dt/2)
		end do
		integer dimension(3) :: f
		total_force(f,gas,n)
		do i = 1,3
			a = -1*r(i)*f(i)/m
		end do
		do i = 1,3
			v(i) = v(i)+.5*a(i)*dt !v(d+dt)
		end do

		do i = 1,3
			gas2(n,1,i)=r(i)
			gas2(n,2,i)=v(i)*m
			gas2(n,3,i)=a(i)*m
		end do
	end do
	gas = gas2

end subroutine verlet_algorithm

subroutine total_force(f,gas,n)!n is the particle we're calculating the force on

	real (kind=ikind), dimension(N,3,3), intent(inout) :: gas
	integer dimension(3), intent(inout) :: f
	integer intent(in) :: n
	

end subroutine gradiant

end program ArgonSim