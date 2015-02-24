module simParam

    ! This module sets all the parameters of the
    ! simulation to be used by other modules

    implicit none

    ! System parameters
    integer :: N
    real(8) :: TEMP
    real(8) :: DENSITY

    ! Simulation parameters
    real(8), parameter :: DT = 4d-3
    integer, parameter :: NUM_STEPS = 10000
    
    ! Temperature renormalization parameters
    integer, parameter :: STEPS_TO_RENORM = 25
    integer, parameter :: NUM_RENORM = 30
    integer, parameter :: FIRST_RENORM  = 100

    ! Size of histogram bin
    real(8), parameter :: DR = 1d-2

    ! Potential cutoff
    real(8), parameter :: CUTOFF = 3.2

    ! Lattice size
    real(8) :: latt_size, box_size
    integer :: N_cube_side

contains
    
    subroutine init_param
    ! Initialize parameters without recurring to file

        ! System parameters
        N = 864
        TEMP = 1
        DENSITY = 0.8

        ! Lattice size
        latt_size = (4d0/DENSITY)**(1d0/3)
        N_cube_side = nint((N/4)**(1d0/3))
        box_size = latt_size*N_cube_side

    end subroutine init_param


  subroutine readSimParam(filename)
  ! Read parameters from a file

    implicit none

    character(len=*), intent(in) :: filename
    character(len=80)  :: token
    integer            :: readStat

    open(unit=20, file=filename, action="read", iostat=readStat)

    do 
      read(20, *, iostat=readStat) token; backspace(20)
      if(readStat .lt. 0) exit

      select case(trim(token))
        case("N")
          read(20, *, iostat=readStat) token, N

        case("TEMP")
          read(20, *, iostat=readStat) token, TEMP

        case("DENSITY")
          read(20, *, iostat=readStat) token, DENSITY

        case default
          read(20, *) !force I/O to advance
      end select

    end do

    ! Print read parameters to screen
    print *, "N: ", N
    print *, "TEMP: ", TEMP
    print *, "DENSITY: ", DENSITY
    print *, "DT: ", DT
    print *, "NUM_STEPS: ", NUM_STEPS
    print *, "STEPS_TO_RENORM: ", STEPS_TO_RENORM
    print *, "NUM_RENORM: ", NUM_RENORM
    print *, "FIRST_RENORM: ", FIRST_RENORM
    print *, "CUTOFF: ", CUTOFF
    print *, "DR: ", DR
    close(20)

    ! Calculate lattice parameters
    latt_size = (4d0/DENSITY)**(1d0/3)
    N_cube_side = nint((N/4)**(1d0/3))
    box_size = latt_size*N_cube_side
    
  end subroutine readSimParam


end module