module verlet

implicit none
public verlet_algorithm, total_force

contains
    subroutine verlet_algorithm(EPS, SIGMA, MASS, l, dt, r, v, a, Ep)
        
      real(8), intent(in) :: dt
      real(8), intent(in) :: l ! size of the box
      real(8), intent(in) :: SIGMA, MASS, EPS
      real (8), intent(inout), dimension(:,:) :: r
      real (8), intent(inout), dimension(:,:) :: v
      real (8), intent(inout), dimension(:,:) :: a
      real (8), intent(out) :: Ep
       
       v = v + 0.5*a*dt
       r = modulo(r + v*dt, l)
       call total_force(EPS, SIGMA, MASS, l, r, a, Ep)
       v = v + 0.5*a*dt

    end subroutine verlet_algorithm

    subroutine total_force(EPS, SIGMA, MASS, l, r, a, Ep) ! n_aux is the particle we're calculating the force on

        real (8), intent(inout), dimension(:,:) :: r
        real (8), intent(inout), dimension(:,:) :: a
        real(8), intent(in) :: SIGMA, MASS, EPS, l
        real(8), intent(out):: Ep
        integer :: N, i, j
        real(8) :: diff_aux(3), dist, fij(3)
        N = size(r,2)
        Ep = 0
        a(:,:) = 0
        do i = 1, N
            do j = (i + 1), N
                diff_aux = r(:, j) - r(:, i)
                diff_aux = diff_aux - nint(diff_aux/l)*l
                dist = sqrt(dot_product(diff_aux, diff_aux))
                fij = (-4d0)* EPS * ((12d0) * (SIGMA/dist) ** (11d0) - &
                 (6d0) * (SIGMA/dist) ** (5d0)) * (diff_aux) / (dist)
                a(:, i) = a(:, i) + fij/MASS
                a(:, j) = a(:, j) - fij/MASS
                Ep = Ep + 4d0*EPS*((SIGMA/dist)**(12d0) - (SIGMA/dist)**(6d0))
            end do
        end do

    end subroutine total_force

end module