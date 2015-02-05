module verlet

implicit none
private total_force
public verlet_algorithm

contains
    subroutine verlet_algorithm(EPS, SIGMA, MASS, l, dt, r, v, a)
        
      real(8), intent(in) :: dt
      real(8), intent(in) :: l ! size of the box
      real(8), intent(in) :: SIGMA, MASS, EPS
      real (8), intent(inout), dimension(:,:) :: r
      real (8), intent(inout), dimension(:,:) :: v
      real (8), intent(inout), dimension(:,:) :: a

       v = v + 0.5*a*dt
       r = modulo(r + v*dt, l)
       call total_force(EPS, SIGMA, MASS, l, r, a)
       v = v + 0.5*a*dt

    end subroutine verlet_algorithm

    subroutine total_force(EPS, SIGMA, MASS, l, r, a) ! n_aux is the particle we're calculating the force on

        real (8), intent(inout), dimension(:,:) :: r
        real (8), intent(inout), dimension(:,:) :: a
        real (8), allocatable, dimension(:,:,:) :: f
        real(8), intent(in) :: SIGMA, MASS, EPS, l
        integer :: N, i, j
        real(8) :: diff_aux(3), dist
        N = size(r,1)
        allocate(f(N,N,3))

        f(:,:,:) = 0
        do i = 1, N
            do j = (i + 1), N
                diff_aux = r(j,:) - r(i, :)
                diff_aux = diff_aux - nint(diff_aux/l)*l
                dist = sqrt(dot_product(diff_aux, diff_aux))
                f(i,j,:) = (-4d0)* EPS * ((12d0) * (SIGMA/dist) ** (11d0) - &
                 (6d0) * (SIGMA/dist) ** (5d0)) * (diff_aux) / (dist)
                f(j,i,:) = -f(i,j,:)
            end do
        end do

        do i = 1, N
            a(i,1) = sum(f(i, :, 1)) / MASS
            a(i,2) = sum(f(i, :, 2)) / MASS
            a(i,3) = sum(f(i, :, 3)) / MASS
        end do
    end subroutine total_force

    subroutine energy(r, v, SIGMA, MASS, EPS, l, Ek, Ep, mom)

      real (8), intent(in), dimension(:,:) :: r
      real (8), intent(in), dimension(:,:) :: v
      real(8), intent(in) :: SIGMA, MASS, EPS, l
      real(8), intent(out) :: Ek, Ep
      real(8), intent(out) :: mom(3)
      integer :: n, i, j
      real(8) :: diff_aux(3), dist


      n = size(v, 1)

      Ek = 0
      Ep = 0
      do i = 1, n
        Ek = Ek + 0.5d0*MASS*dot_product(v(i,:), v(i,:))
      end do
      open(unit = 1000, file = "dist.dat") 
      do i = 1, n
        do j = i+1, n
          diff_aux = r(i,:) - r(j, :)
          diff_aux = diff_aux - nint(diff_aux/l)*l
          dist = sqrt(dot_product(diff_aux, diff_aux))
          write(1000, *) dist
          Ep = Ep + 4d0*EPS*((SIGMA/dist)**(12d0) - (SIGMA/dist)**(6d0))
        end do
      end do
      close(1000)
      do i = 1, 3
        mom(i) = MASS * sum(v(:,i))
      end do
    end subroutine energy
end module