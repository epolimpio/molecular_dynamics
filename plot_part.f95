MODULE Plot_Part
  implicit none

contains

  subroutine Init_Plot

    use simParam, only: box_size

    call InitPlot('lightblue', 700, 700, 'out', 1)
    call Framing(0._8, 0._8, box_size, box_size)
    call PutStopButton()

  end subroutine Init_Plot

  subroutine Plot_Particles(Pos, N)
      integer :: N, i
      real(8) :: Pos(3,N)
      do i = 1, N
        call setpoint(Pos(1,i), Pos(2,i))
      end do
  end subroutine Plot_Particles

end module plot_part
  
