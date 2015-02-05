program ArgonSim

    use initial_conditions
    use verlet
    !use sub_regions
    implicit none

    !integer, parameter  :: ikind=selected_real_kind(p=15)
    
    ! System parameters
    integer, parameter :: N = 108
    real(8), parameter :: SIGMA = 1
    real(8), parameter :: EPS = 1
    real(8), parameter :: MASS = 1
    real(8), parameter :: T = 1
    !real(8), parameter :: VOL = 3
    real(8), parameter :: DT = 1d-3

    ! Initialize array structure
    real (8), dimension(N,3) :: r
    real (8), dimension(N,3) :: v
    real (8), dimension(N,3) :: a
    real(8) :: latt_size, cube_side, l, Ek, Ep
    real(8) :: mom(3)
    integer(8) :: i, k

    a(:,:) = 0

    latt_size = 2d0**(2d0/3) * SIGMA
    cube_side = nint((N/4)**(1d0/3))

    l = latt_size*cube_side

    ! Initialize random distribution(s) & lattices

    call max_boltz(MASS, T, v)
    call fcc_lattice(SIGMA, r)
    call energy(r, v, SIGMA, MASS, EPS, l, Ek, Ep, mom)

    print *, "Initial Momenta: ", mom
    print *, "Initial Ek: ", Ek
    print *, "Initial Ep: ", Ep
    print *, "Initial E: ", Ek+Ep

    do i = 1, 100
        call verlet_algorithm(EPS, SIGMA, MASS, l, DT, r, v, a)
        call energy(r, v, SIGMA, MASS, EPS, l, Ek, Ep, mom)
        print *, "Ek: ", Ek, "Ep: ", Ep, "E: ", Ek+Ep
     
    end do

    !do currParticle = 1,N
    !    ! TODO Remove momentum check if it works
    !    ! Make sure Momentums sum to 0 at end of loop
    !    if (currParticle == N) then
    !        do currAxis = 1,3
    !            if (totalMomentum(currAxis) > 0.00001) then
    !            print *, "Momentum NOT conserved!"
    !                ! gas(currentParticle,Momentums,currAxis) = -totalMomentum(currAxis)
    !            end if
    !        end do
    !    end if
    !    ! Move current particle based on current momentum
    !    do currAxis = 1,3
    !        !gas(currParticle,Positions,currAxis) = gas(currParticle,Positions,currAxis) + (gas(currParticle,Momentums,currAxis)/mass)*dt
    !    end do
    !    ! Increment time
    !end do
    

end program ArgonSim

