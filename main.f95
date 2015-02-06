program ArgonSim

    use initial_conditions
    use verlet
    use physical_quantities
    use plot_part
    !use sub_regions
    implicit none
    
    ! System parameters
    integer, parameter :: N = 864
    real(8), parameter :: SIGMA = 1
    real(8), parameter :: EPS = 1
    real(8), parameter :: MASS = 1
    real(8), parameter :: T = 1.3
    real(8), parameter :: DENSITY = 0.8

    ! Simulation parameters
    real(8), parameter :: DT = 4d-3
    integer, parameter :: NUM_STEPS = 2000
    
    ! Temperature renormalization parameters
    integer, parameter :: STEPS_TO_RENORM = 30
    real(8), parameter :: T_DRIFT = 5d-3

    ! Initialize array structure
    real (8), dimension(3, N) :: r
    real (8), dimension(3, N) :: v
    real (8), dimension(3, N) :: a
    real(8) :: latt_size, cube_side, l, Ek, Ep, temp_out
    real(8) :: mom(3)
    integer(8) :: i, frame, k
    character(5) :: strframe
    logical :: renorm
    
    ! Counter frame used to generate positions data
    frame = 0

    ! Start with renormalization
    renorm = .TRUE.

    ! Determine the volume of the system and its size
    latt_size = (4d0/DENSITY)**(1d0/3) * SIGMA
    cube_side = nint((N/4)**(1d0/3))
    l = latt_size*cube_side

    ! Initialize random distributions & lattices
    call max_boltz(MASS, T, v)
    call fcc_lattice(latt_size, r)
    
    ! Initialize graphics
    call InitPlot('lightblue', 700, 700, 'out', 1)
    call Framing(0._8, 0._8, l, l)
    call PutStopButton()

    ! Calculate initial energies and initialize accelerations
    call kinetic_energy(v, MASS, Ek, mom, temp_out)
    call total_force(EPS, SIGMA, MASS, l, r, a, Ep)

    open(unit = 1001, file = "data/energy.dat")
    open(unit = 1002, file = "data/temperature.dat")

    ! Print initial data
    print *, "Initial Momenta: ", mom
    print *, "Initial Ek: ", Ek
    print *, "Initial Ep: ", Ep
    print *, "Initial E: ", Ek+Ep
    do i = 1, NUM_STEPS
        ! Runs the verlet time evolution for each time step
        call verlet_algorithm(EPS, SIGMA, MASS, l, DT, r, v, a, Ep)
        call kinetic_energy(v, MASS, Ek, mom, temp_out)
        write(1001, *) Ek, Ep, Ek+Ep, temp_out
        call plot_particles(r, N)
        ! Check for renormalization and renormalize
        if (mod(i, STEPS_TO_RENORM) == 0 .AND. renorm) then
            ! Check if temperature has drift much
            if (abs(T - temp_out) < T_DRIFT) then
                print *, "Initial Momenta: ", mom
                renorm = .FALSE.
            else
                v = sqrt(1d0*T/temp_out)*v
            end if
        end if
        ! Monitor the results printing positions and energies to files
        if (mod(i,100) == 0) then
            frame = frame + 1
            print *, frame, "Ek: ", Ek, "Ep: ", Ep, "E: ", Ek+Ep, "T: ", temp_out
            write(strframe,'(I5)') frame
            open(unit = 1000, file = "data/position"//strframe//".dat")
            do k = 1, N
              write(1000, *) r(:,k)
            end do
           close(1000)
        end if

    end do

    call EndPlot()
    close(1001)

end program ArgonSim

