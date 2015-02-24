module simulation

use simParam
use plot_part
use physical_quantities
use output_file

implicit none

real (8), dimension(:, :), allocatable :: r
real (8), dimension(:, :), allocatable :: v
real (8), dimension(:, :), allocatable :: a
real (8) :: Ep, Ek, temp_out, sum_rF

contains

  subroutine init_all(init_parameters)
  ! This subroutine initializes all varibles
  ! It reads the parameters from a file
  ! and initialize positions and velocities

    use initial_conditions

    logical, intent(in) :: init_parameters

    ! Read parameters from a file
    if (init_parameters) then
      call readSimParam("./input.ini")
    end if

    ! Initialize correlation function histogram
    ! and calculation of heat capacity
    call init_histogram

    ! Allocate the particle dynamics vectors
    allocate(r(3,N))
    allocate(v(3,N))
    allocate(a(3,N))

    ! Initialize positions and velocities
    ! and calculate initial energies
    call initial_values(r, v)
    call total_force(.FALSE.)
    call kinetic_energy
    print *, "Initial Momenta: ", sum(v,2)
    print *, "Initial Ek: ", Ek
    print *, "Initial Ep: ", Ep
    print *, "Initial E: ", Ek+Ep

  end subroutine init_all

  subroutine stabilize_temp
  ! This subroutine runs the dynamics with the termostat
  ! until realeased with stable temperature.
  ! This is done by multiplying by a factor.

    integer(8) :: i, j
    
    ! Runs the dynamics until the step of first renormalization
    do i = 1, FIRST_RENORM
      call update_verlet(.FALSE.)
    end do
    call kinetic_energy
      
    ! Print Results
    print *, "1st Renorm... ", "Ek: ", Ek, "Ep: ", Ep, "E: ", Ek+Ep, "T: ", temp_out

    ! First Renormalization
    v = sqrt(1d0*TEMP/temp_out)*v    

    ! Runs the dynamics, renormalizing every STEPS_TO_RENORM.
    ! This is done NUM_RENORM times
    do i = 1, NUM_RENORM
      do j = 1, STEPS_TO_RENORM
        call update_verlet(.FALSE.)
      end do
      call kinetic_energy
      
      ! Print Results
      print *, "Renorm... ", "Ek: ", Ek, "Ep: ", Ep, "E: ", Ek+Ep, "T: ", temp_out
      
      ! Renormalize
      v = sqrt(1d0*TEMP/temp_out)*v
    end do

  end subroutine

  subroutine simulate_dynamics
  ! This subroutine simulate the dynamics according to
  ! the number of steps (NUM_STEPS) provided
    integer(8) :: i

    ! Initializes plotting to see the particles movement
    ! call Init_Plot

    do i = 1, NUM_STEPS
      
      ! The .TRUE. in update_verlet means that distances
      ! must be added to correlation histogram 
      call update_verlet(.TRUE.)

      ! Plot particle movement and calc energy
      ! call plot_particles(r, N)
      call kinetic_energy

      ! Save parameters to calc physical quantities
      call save_kinetic(i, Ek)
      call save_forces(i, sum_rF)
      call save_velocities(i, v)

      ! Write energies to a file
      call write_energies(Ek, Ep, temp_out)

      ! Print results on screen every 100 steps
      if (mod(i,100) == 0) then
        print *, i, "Ek: ", Ek, "Ep: ", Ep, "E: ", Ek+Ep, "T: ", temp_out
      end if
    end do

    ! Normalize distances histogram to calculate
    ! the correlation function
    call normalize_histogram
    ! Calculate heat capacity
    call calc_heat_capacity
    call calc_pressure
    call calc_veloc_corr

    print *, "Cv: ", Cv, "P: ", P

    ! Write correlation function to a file
    call write_correlation
    call write_vel_corr

    print *, "END!"

    ! Close plotting
    ! call EndPlot()

  end subroutine simulate_dynamics

  subroutine update_verlet(calc_correl)
  ! Update positions according to Verlet algorithm.
  ! calc_correl indicates if distances must be added
  ! to the correlation function in force calculation.

    logical, intent(in) :: calc_correl
      
    v = v + 0.5*a*DT
    ! The modulo deals with the boudary conditions
    r = modulo(r + v*DT, box_size)
    call total_force(calc_correl)
    v = v + 0.5*a*DT

  end subroutine update_verlet

  subroutine total_force(calc_correl)
  ! Calculate the total forces of the particles
  ! It also calculates the potential energy and 
  ! adds distances to the correlation function
  ! according to the parameter calc_correl

    logical, intent(in) :: calc_correl
    integer :: i, j
    real(8) :: diff_aux(3), dist2, fij, dist
    
    ! Initialize variables
    Ep = 0
    a(:,:) = 0
    sum_rF = 0

    ! Calculate for all pairs
    do i = 1, N
        do j = (i + 1), N
            ! Vector difference and distance squared
            diff_aux = r(:, i) - r(:, j)
            ! This takes into account the nieghbour cells
            diff_aux = diff_aux - nint(diff_aux/box_size)*box_size
            dist2 = dot_product(diff_aux, diff_aux)
            dist = sqrt(dist2)
            
            ! Calculate if the distance is below cutoff
            if (dist2 .LE. CUTOFF**2d0) then
              fij = (4d0)* ((12d0) * (1d0/dist2) ** (7d0) - &
               (6d0) * (1d0/dist2) ** (4d0))
              a(:, i) = a(:, i) + diff_aux*fij
              a(:, j) = a(:, j) - diff_aux*fij
              sum_rF = sum_rF + dist2*fij
            end if

            Ep = Ep + 4d0*((1d0/dist2)**(6d0) - (1d0/dist2)**(3d0))

            ! Place distances on correlation func histogram
            if (calc_correl) then
              call place_on_histogram(dist)
            end if 
        end do
    end do

  end subroutine total_force

  subroutine kinetic_energy
  ! Calculate kinetic energy and the corresponding
  ! temperature, using Boltzmann const = 1

    integer :: i
    Ek = 0
    do i = 1, N
      Ek = Ek + 0.5d0*dot_product(v(:, i), v(:, i))
    end do
    temp_out = 2d0*Ek/(3d0*(n-1))

  end subroutine kinetic_energy 

end module