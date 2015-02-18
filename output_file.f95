module output_file
    
    use simParam, only: FNAME, TEMP, N, DENSITY
    implicit none

    integer, parameter :: fEnergy = 1001
    integer, parameter :: fCorrel = 1002
    integer, parameter :: fHeatCap = 1003
    !integer, parameter :: fPressure = 1004

contains

    subroutine init_files

        open(unit = fEnergy, file = "./data/energy" // FNAME // ".dat")
        open(unit = fCorrel, file = "./data/correlation" // FNAME // ".dat")
        open(unit = fHeatCap, file = "./data/heat_cap.dat", access = 'APPEND')

        call insert_header(fEnergy)
        call insert_header(fCorrel)


    end subroutine init_files

    subroutine insert_header(file)

        integer, intent(in) :: file

        write(file, *) "TEMP->", TEMP
        write(file, *) "DENSITY->", DENSITY
        write(file, *) "N->", N

    end subroutine

    subroutine close_files

        close(fEnergy)
        close(fCorrel)
        close(fHeatCap)

    end subroutine close_files


end module