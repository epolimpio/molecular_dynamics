program md

  use initial_conditions

  implicit none
  integer, parameter :: N = 256
  real(8), parameter :: T = 1
  real(8) :: a(N,3), r(N,3), v(N,3)
  integer :: i

  open(unit = 1000, file = "momentum.dat")
  open(unit = 1001, file = "position.dat")
  call max_boltz(T, v)
  call fcc_lattice(1d0, r)
  do i = 1, N
    write(1000, *) v(i, 1)
    write(1001, *) r(i, :)
  end do
  
  close(1000)

end program
  
  
  

	
