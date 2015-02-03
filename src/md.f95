program md

  use initial_conditions

  implicit none

  integer, parameter :: N = 2048
  real(8), parameter :: T = 1
  real(8) :: momenta(3, N)
  integer :: i

  open(unit = 1000, file = "momentum.dat")
  call max_boltz(N, T, momenta)
  do i = 1, N
    write(1000, *) momenta(1,i)
  end do
  
  close(1000)

end program
  
  
  

	
