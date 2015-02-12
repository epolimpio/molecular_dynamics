module verlet

use simParam

implicit none

real (8), dimension(:, :), allocatable :: r
real (8), dimension(:, :), allocatable :: v
real (8), dimension(:, :), allocatable :: a
real (8) :: Ep, Ek, temp_out

contains

  subroutine init_all

    use initial_conditions

    call init_param
    allocate(r(3,N))
    allocate(v(3,N))
    allocate(a(3,N))
    call initial_values(r, v)
    call total_force
    call kinetic_energy
    print *, "Initial Momenta: ", sum(v,2)
    print *, "Initial Ek: ", Ek
    print *, "Initial Ep: ", Ep
    print *, "Initial E: ", Ek+Ep

  end subroutine init_all

  subroutine stabilize_temp

    integer(8) :: i, j
    do i = 1, NUM_RENORM
      do j = 1, STEPS_TO_RENORM
        call update_verlet
      end do
      call kinetic_energy
      v = sqrt(1d0*TEMP/temp_out)*v
    end do

  end subroutine

  subroutine simulate_dynamics

    integer(8) :: i

    do i = 1, NUM_STEPS
      
      call update_verlet
      call kinetic_energy
      
      if (mod(i,10) == 0) then
        print *, i, "Ek: ", Ek, "Ep: ", Ep, "E: ", Ek+Ep, "T: ", temp_out
      end if
    end do

  end subroutine simulate_dynamics

  subroutine update_verlet
      
    v = v + 0.5*a*DT
    r = modulo(r + v*DT, box_size)
    call total_force
    v = v + 0.5*a*DT

  end subroutine update_verlet

  subroutine total_force

    integer :: i, j
    real(8) :: diff_aux(3), dist, fij(3)
        
    Ep = 0
    a(:,:) = 0
    do i = 1, N
        do j = (i + 1), N
            diff_aux = r(:, j) - r(:, i)
            diff_aux = diff_aux - nint(diff_aux/box_size)*box_size
            dist = sqrt(dot_product(diff_aux, diff_aux))
            fij = (-4d0)* ((12d0) * (1d0/dist) ** (11d0) - &
             (6d0) * (1d0/dist) ** (5d0)) * (diff_aux) / (dist)
            a(:, i) = a(:, i) + fij
            a(:, j) = a(:, j) - fij
            Ep = Ep + 4d0*((1d0/dist)**(12d0) - (1d0/dist)**(6d0))
        end do
    end do

  end subroutine total_force

  subroutine kinetic_energy

    integer :: i
    Ek = 0
    do i = 1, N
      Ek = Ek + 0.5d0*dot_product(v(:, i), v(:, i))
    end do
    temp_out = 2d0*Ek/(3d0*(n-1))

  end subroutine kinetic_energy 

end module