program ArgonSim

    use initial_conditions
    use verlet
    !use sub_regions
    implicit none

    !integer, parameter  :: ikind=selected_real_kind(p=15)
    
    ! System parameters
    integer, parameter :: N = 864
    real(8), parameter :: SIGMA = 1
    real(8), parameter :: EPS = 1
    real(8), parameter :: MASS = 1
    real(8), parameter :: T = 1
    real(8), parameter :: DENSITY = 0.88
    real(8), parameter :: DT = 1d-3

    ! Initialize array structure
    real (8), dimension(3, N) :: r
    real (8), dimension(3, N) :: v
    real (8), dimension(3, N) :: a
    real(8) :: latt_size, cube_side, l, Ek, Ep
    real(8) :: mom(3)
    integer(8) :: i, frame, k
    character(5) :: strframe
    
    ! Counter frame used to generate positions data
    frame = 0

    ! Determine the volume of the system and its size
    latt_size = 2d0**(2d0/3) * SIGMA
    cube_side = nint((N/4)**(1d0/3))
    l = 1.0d0*latt_size*cube_side

    ! Initialize random distributions & lattices
    call max_boltz(MASS, T, v)
    call fcc_lattice(SIGMA, r)

    ! Calculate initial energies and initialize accelerations
    call kinetic_energy(v, MASS, Ek, mom)
    call total_force(EPS, SIGMA, MASS, l, r, a, Ep)

    open(unit = 1001, file = "data/energy.dat")

    ! Print initial data
    print *, "Initial Momenta: ", mom
    print *, "Initial Ek: ", Ek
    print *, "Initial Ep: ", Ep
    print *, "Initial E: ", Ek+Ep
    do i = 1, 10000
        ! Runs the verlet time evolution for each time step
        call verlet_algorithm(EPS, SIGMA, MASS, l, DT, r, v, a, Ep)
        call kinetic_energy(v, MASS, Ek, mom)
        write(1001, *) Ek, Ep, Ek+Ep
        
        ! Monitor the results printing positions and energies to files
        if (mod(i,100) == 0) then
            frame = frame + 1
            print *, frame, "Ek: ", Ek, "Ep: ", Ep, "E: ", Ek+Ep
            write(strframe,'(I5)') frame
            open(unit = 1000, file = "data/position"//strframe//".dat")
            do k = 1, N
              write(1000, *) r(k,:)
            end do
           close(1000)
        end if

    end do

    close(1001)

end program ArgonSim

