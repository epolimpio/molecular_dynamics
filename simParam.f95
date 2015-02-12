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
    integer :: FIRST_RENORM

    ! Size of histogram bin
    real(8) :: DR

    ! Lattice size
    real(8) :: latt_size, N_cube_side, box_size

contains
    
    subroutine init_param

        ! System parameters
        N = 864
        TEMP = 1
        DENSITY = 0.8

        ! Simulation parameters
        DT = 4d-3
        NUM_STEPS = 1500

        ! Temperature renormalization parameters
        FIRST_RENORM = 200
        STEPS_TO_RENORM = 2
        NUM_RENORM = 150

        ! Potential cutoff
        CUTOFF = 4 ! 4*sigma

        ! Histogram bin size
        DR = 1d-3 ! 1/100 of sigma

        ! Lattice size
        latt_size = (4d0/DENSITY)**(1d0/3)
        N_cube_side = nint((N/4)**(1d0/3))
        box_size = latt_size*N_cube_side

    end subroutine init_param

end module