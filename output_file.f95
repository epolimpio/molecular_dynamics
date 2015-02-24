module output_file
    
    use simParam
    use physical_quantities, only: correl_histogram, &
                            vel_corr, vel_corr_err, &
                            mod_vel_corr, mod_vel_corr_err
    implicit none

    integer, parameter :: fEnergy = 1001
    integer, parameter :: fCorrel = 1002
    integer, parameter :: fVelCorr = 1003
    integer, parameter :: fPressure = 1004

contains

    subroutine init_files(fname)

        character(len=*), intent(in) :: fname

        open(unit = fEnergy, file = "./data/energy" // fname // ".dat")
        open(unit = fCorrel, file = "./data/correlation" // fname // ".dat")
        open(unit = fPressure, file = "./data/pressure" // fname // ".dat")
        open(unit = fVelCorr, file = "./data/vel_corr" // fname // ".dat")

        call insert_header(fEnergy)
        call insert_header(fCorrel)
        call insert_header(fPressure)
        call insert_header(fVelCorr)

    end subroutine init_files

    subroutine insert_header(file)

        integer, intent(in) :: file

        write(file, *) "TEMP->", TEMP
        write(file, *) "DENSITY->", DENSITY
        write(file, *) "N->", N
        write(file, *) "CUTOFF->", CUTOFF
        write(file, *) "NUM_STEPS->", NUM_STEPS

    end subroutine

    subroutine write_energies(Ek, Ep, T)

        real(8), intent(in) :: Ek, Ep, T

        write(fEnergy, *) Ek, Ep, Ek+Ep, T

    end subroutine write_energies

    subroutine write_correlation

        integer :: i

        do i = 1, size(correl_histogram)
          write(fCorrel, *) correl_histogram(i)
        end do

    end subroutine write_correlation

    subroutine write_vel_corr

        integer :: i

        do i = 1, size(vel_corr)
          write(fVelCorr, *) vel_corr(i), vel_corr_err(i), mod_vel_corr(i), mod_vel_corr_err(i)
        end do

    end subroutine write_vel_corr

    subroutine close_files

        close(fEnergy)
        close(fCorrel)
        close(fPressure)
        close(fVelCorr)

    end subroutine close_files


end module