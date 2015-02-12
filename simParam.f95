module simParam

    ! This module sets all the parameters of the
    ! simulation to be used by other modules

    implicit none

    ! System parameters
    integer :: N
    real(8) :: TEMP
    real(8) :: DENSITY

    ! Simulation parameters
    real(8) :: DT
    integer :: NUM_STEPS
    
    ! Temperature renormalization parameters
    integer :: STEPS_TO_RENORM
    integer :: NUM_RENORM

    ! Lattice size
    real(8) :: latt_size, N_cube_side, box_size

contains
    
    subroutine init_param

        ! System parameters
        N = 864
        TEMP = 1.3
        DENSITY = 1

        ! Simulation parameters
        DT = 4d-3
        NUM_STEPS = 2000

        ! Temperature renormalization parameters
        STEPS_TO_RENORM = 30
        NUM_RENORM = 20

        ! Lattice size
        latt_size = (4d0/DENSITY)**(1d0/3)
        N_cube_side = nint((N/4)**(1d0/3))
        box_size = latt_size*N_cube_side

    end subroutine init_param

end module