module verlet

use simParam
use plot_part
use physical_quantities

implicit none

real (8), dimension(:, :), allocatable :: r
real (8), dimension(:, :), allocatable :: v
real (8), dimension(:, :), allocatable :: a
real (8) :: Ep, Ek, temp_out

contains

  subroutine init_all

    use initial_conditions

    call init_param
    call init_histogram
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

    call Init_Plot

  end subroutine init_all

  subroutine stabilize_temp

    integer(8) :: i, j
    
    do i = 1, FIRST_RENORM
      call update_verlet
    end do
    call kinetic_energy
      
    ! Print Results
    print *, "1st Renorm... ", "Ek: ", Ek, "Ep: ", Ep, "E: ", Ek+Ep, "T: ", temp_out

    ! Renormalize
    v = sqrt(1d0*TEMP/temp_out)*v    

    do i = 1, NUM_RENORM
      do j = 1, STEPS_TO_RENORM
        call update_verlet
      end do
      call kinetic_energy
      
      ! Print Results
      print *, "Renorm... ", "Ek: ", Ek, "Ep: ", Ep, "E: ", Ek+Ep, "T: ", temp_out
      
      ! Renormalize
      v = sqrt(1d0*TEMP/temp_out)*v
    end do

  end subroutine

  subroutine simulate_dynamics

    integer(8) :: i

    open(unit = 1001, file = "data/energy.dat")
    open(unit = 1002, file = "data/correlation.dat")

    do i = 1, NUM_STEPS
      
      call update_verlet
      call kinetic_energy

      write(1001, *) Ek, Ep, Ek+Ep, temp_out

      ! Print results every 100 steps
      if (mod(i,100) == 0) then
        print *, i, "Ek: ", Ek, "Ep: ", Ep, "E: ", Ek+Ep, "T: ", temp_out
      end if
    end do

    call normalize_histogram

    do i = 1, size(correl_histogram)
      write(1002, *) correl_histogram(i)
    end do

    close(1001)
    close(1002)

    print *, "END!"

    call EndPlot()

  end subroutine simulate_dynamics

  subroutine update_verlet
      
    v = v + 0.5*a*DT
    r = modulo(r + v*DT, box_size)
    call plot_particles(r, N)
    call total_force
    v = v + 0.5*a*DT

  end subroutine update_verlet

  subroutine total_force

    integer :: i, j
    real(8) :: diff_aux(3), dist2, fij
        
    Ep = 0
    a(:,:) = 0
    do i = 1, N
        do j = (i + 1), N
            diff_aux = r(:, j) - r(:, i)
            diff_aux = diff_aux - nint(diff_aux/box_size)*box_size
            dist2 = dot_product(diff_aux, diff_aux)

            if (dist2 .LE. CUTOFF**2d0) then
              fij = (-4d0)* ((12d0) * (1d0/dist2) ** (7d0) - &
               (6d0) * (1d0/dist2) ** (4d0))
              a(:, i) = a(:, i) + diff_aux*fij
              a(:, j) = a(:, j) - diff_aux*fij
              Ep = Ep + 4d0*((1d0/dist2)**(6d0) - (1d0/dist2)**(3d0))
            end if
            
            call place_on_histogram(sqrt(dist2))
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