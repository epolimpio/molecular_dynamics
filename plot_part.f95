MODULE Plot_Part
  implicit none

contains

  subroutine Plot_Particles(Pos, N)
  integer :: N, i
  real(8) :: Pos(3,N)
  do i = 1, N
    call setpoint(Pos(1,i), Pos(2,i))
  end do
  end subroutine Plot_Particles
end module plot_part
  
