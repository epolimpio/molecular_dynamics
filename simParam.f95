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

    ! Potential cutoff
    real(8) :: CUTOFF

    character(2) :: FNAME

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

        ! Simulation parameters
        DT = 4d-3
        NUM_STEPS = 1500

        ! Temperature renormalization parameters
        FIRST_RENORM = 100
        STEPS_TO_RENORM = 40
        NUM_RENORM = 10

        ! Potential cutoff
        CUTOFF = 4 ! 4*sigma

        ! Histogram bin size
        DR = 1d-3 ! 1/100 of sigma

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

        case("DT")
          read(20, *, iostat=readStat) token, DT

        case("NUM_STEPS")
          read(20, *, iostat=readStat) token, NUM_STEPS

        case("FIRST_RENORM")
          read(20, *, iostat=readStat) token, FIRST_RENORM

        case("STEPS_TO_RENORM")
          read(20, *, iostat=readStat) token, STEPS_TO_RENORM

        case("NUM_RENORM")
          read(20, *, iostat=readStat) token, NUM_RENORM

        case("CUTOFF")
          read(20, *, iostat=readStat) token, CUTOFF

        case("DR")
          read(20, *, iostat=readStat) token, DR

        case("FNAME")
          read(20, *, iostat=readStat) token, FNAME

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
    print *, "FNAME: ", "./data/teste" // FNAME // ".txt"
    close(20)

    ! Calculate lattice parameters
    latt_size = (4d0/DENSITY)**(1d0/3)
    N_cube_side = nint((N/4)**(1d0/3))
    box_size = latt_size*N_cube_side
    
  end subroutine readSimParam

end module