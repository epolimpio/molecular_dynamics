module physical_quantities

  implicit none
  private
  public kinetic_energy

contains

    subroutine kinetic_energy(v, MASS, Ek, mom, temp_out)

      real (8), intent(in), dimension(:,:) :: v
      real(8), intent(in) :: MASS
      real(8), intent(out) :: Ek, temp_out
      real(8), intent(out) :: mom(3)
      integer :: n, i

      n = size(v, 2)

      Ek = 0
      do i = 1, n
        Ek = Ek + 0.5d0*MASS*dot_product(v(:, i), v(:, i))
      end do
      temp_out = 2d0*Ek/(3d0*(n-1))
      do i = 1, 3
        mom(i) = MASS * sum(v(i, :))
      end do
    end subroutine kinetic_energy

end module