module physical_quantities
    
    use simParam, only: box_size, DR, N, NUM_STEPS

    implicit none
    real(8), parameter :: PI = 4d0*atan(1d0)
    real(8), allocatable :: correl_histogram(:)
    real(8), allocatable :: E_kin(:)
    real(8) :: Cv
    integer :: num_bins

contains

    subroutine init_histogram

        ! We get the number of bins
        ! +1 to get the tail (distances above box_size/2)
        ! The maximum distance in a specific direction is
        ! half the box size, we choose sqrt times that as max. bin
        num_bins = floor(box_size/2d0/DR) + 1
        allocate(correl_histogram(num_bins))
        correl_histogram(:) = 0

    end subroutine init_histogram

    subroutine place_on_histogram(dist)

        real(8), intent(in) :: dist
        integer :: index

        index = floor(dist/DR) + 1
        if (index .LE. num_bins) then
            correl_histogram(index) = correl_histogram(index) + 1
        end if

    end subroutine place_on_histogram

    subroutine normalize_histogram

        integer :: i
        real(8) :: radius

        ! We first divide by the number of runs
        correl_histogram(:) = correl_histogram(:)/NUM_STEPS

        ! We then normalize by the density
        correl_histogram(:) = 2d0*(box_size**3d0)*correl_histogram(:)/(N*(N-1))

        do i = 1, num_bins
            ! We use the radius corresponding to the middle of the bin
            radius = i*DR - DR/2
            correl_histogram(i) = correl_histogram(i)/(DR*4*PI*radius**2d0)
        end do

    end subroutine normalize_histogram

    subroutine save_kinetic(i, Ek)

        real(8), intent(in) :: Ek
        integer(8), intent(in) :: i

        ! Allocate the vector if not done
        if (.NOT. allocated(E_kin)) allocate(E_kin(NUM_STEPS))

        E_kin(i) = Ek

    end subroutine save_kinetic

    subroutine calc_heat_capacity

        real(8) :: mean_Ek
        real(8) :: dev_Ek_sq

        mean_Ek = sum(E_kin)/NUM_STEPS
        dev_Ek_sq = sum((E_kin - mean_Ek)**2d0)/NUM_STEPS

        Cv = 3d0*N/2 / (1 - (3d0*N*dev_Ek_sq/2/(mean_Ek**2)))

    end subroutine calc_heat_capacity

end module