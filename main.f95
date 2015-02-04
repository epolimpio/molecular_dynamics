program ArgonSim
    
    ! TODO ----------------------------------
    !       -Compartmentalize code
    !       -Find and utilize visualization plugin
    ! ---------------------------------------

    use initial_conditions
    implicit none

    ! Initialize misc variables
    integer :: forces,momentums,positions,x,y,z
    positions = 1 ! To improve code usability; e.g. gas(1,forces,x)...
    forces = 2    ! ...or gas(242,momentums,z)
    momentums = 3
    x = 1
    y = 2
    z = 3

    integer :: N,currParticle,currKinem,currAxis,compareParticle
    real, dimension(3) :: totalMomentum
    N = 1024

    ! Initialize array structure
    real (kind=ikind), dimension(N,3,3) :: gas
    integer, allocable :: relatedParticles(N,:)
    real :: compareDistance, cutDistance
    cutDistance = 1

    ! Initialize random distribution(s)

    ! Fill array with initial values
        ! Positions based on FCC
        ! Forces calculated based on positions ?
        ! Momentums random normal dist
    do currParticle = 1,N
        do currKinem = 1,3
            do currAxis = 1,3

            end do
        end do
    end do

    ! Loop over gas
    do currParticle = 1,N
        ! TODO Remove momentum check if it works
        ! Make sure momentums sum to 0 at end of loop
        if (currParticle == N) then
            do currAxis = 1,3
                if (totalMomentum(currAxis) > 0.00001) then
                print *, "Momentum NOT conserved!"
                    ! gas(currentParticle,momentums,currAxis) = -totalMomentum(currAxis)
                end if
            end do
        end if
        ! Move current particle based on current momentum
        do currAxis = 1,3
            !gas(currParticle,positions,currAxis) = gas(currParticle,positions,currAxis) + (gas(currParticle,momentums,currAxis)/mass)*dt
        end do
        ! Increment time
    end do

end program ArgonSim

