module physical_quantities
    
    use simParam

    implicit none
    real(8), parameter :: PI = 4d0*atan(1d0)
    real(8), allocatable :: correl_histogram(:)
    real(8), allocatable :: E_kin(:)
    real(8), allocatable :: rF(:)
    real(8), allocatable :: dr_sq(:,:)
    real(8), allocatable :: vel(:,:,:)
    real(8), allocatable :: vel_corr(:), mod_vel_corr(:)
    real(8), allocatable :: vel_corr_err(:), mod_vel_corr_err(:)
    real(8) :: Cv, P
    integer :: num_bins

    integer, parameter :: NUM_DT = 1000

contains

    subroutine init_histogram

        ! We get the number of bins
        ! +1 to get the tail (distances above box_size/2)
        ! The maximum distance in a specific direction is
        ! half the box size, we choose sqrt times that as max. bin
        num_bins = floor(box_size/2d0/DR) + 1
        if (allocated(correl_histogram)) deallocate(correl_histogram)
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

        print *, "Av. Ek", mean_Ek, "+-", sqrt(dev_Ek_sq)

        Cv = (1d0/(1.5d0*N)-dev_Ek_sq/mean_Ek**2d0)**(-1d0)

    end subroutine calc_heat_capacity

    subroutine save_forces(i, sum_rF)

        real(8), intent(in) :: sum_rF
        integer(8), intent(in) :: i

        ! Allocate the vector if not done
        if (.NOT. allocated(rF)) allocate(rF(NUM_STEPS))

        rF(i) = sum_rF

    end subroutine save_forces

    subroutine calc_pressure

        P = 1 + 1d0/(3d0*N*TEMP)*sum(rF)/NUM_STEPS &
            + 16d0*PI*DENSITY/(3d0*TEMP)*(2d0/(3d0*CUTOFF**9) - 1d0/(CUTOFF**3))

    end subroutine calc_pressure

    subroutine save_velocities(i, velocities)

        real(8), intent(in) :: velocities(:, :)
        integer(8), intent(in) :: i

        ! Allocate the vector if not done
        if (.NOT. allocated(vel)) allocate(vel(3, N, NUM_STEPS))

        vel(:, :, i) = velocities

    end subroutine save_velocities

    subroutine calc_veloc_corr

        integer :: tau, i, count, part
        real(8) :: corr, corr_err, mean, mean_mod, mod_corr, mod_corr_err
        real(8) :: mean_v, var_v, v(N, NUM_STEPS)

        ! Allocate vectors
        if (.NOT. allocated(vel_corr)) allocate(vel_corr(0:NUM_DT))
        if (.NOT. allocated(vel_corr_err)) allocate(vel_corr_err(0:NUM_DT))

        if (.NOT. allocated(mod_vel_corr)) allocate(mod_vel_corr(0:NUM_DT))
        if (.NOT. allocated(mod_vel_corr_err)) allocate(mod_vel_corr_err(0:NUM_DT))

        ! Calculate velocity modulus
        do i = 1, NUM_STEPS
            do part = 1, N
                v(part, i) = sqrt(dot_product(vel(:, part, i), vel(:, part, i)))
            end do
        end do

        mean_v = sum(sum(v,1)/N)/NUM_STEPS
        v = v - mean_v
        var_v = sum(sum(v**2d0,1)/N)/NUM_STEPS

        print *, "Mean vel.: ", mean_v, "+-", sqrt(var_v)
        
        do tau = 0, NUM_DT
            i = 1
            count = 0
            corr = 0
            corr_err = 0
            mod_corr = 0
            mod_corr_err = 0
            do while (i+tau .LE. NUM_STEPS)
                mean = 0
                mean_mod = 0
                count = count + 1
                do part  = 1, N
                    mean = mean + dot_product(vel(:, part, i), vel(:, part, i+tau))
                    mean_mod = mean_mod + v(part,i)*v(part, i+tau)
                end do
                mean = 1d0*mean/N
                mean_mod = 1d0*mean_mod/N
                corr = corr + mean
                corr_err = corr_err + mean**2d0
                mod_corr = mod_corr + mean_mod
                mod_corr_err = mod_corr_err  + mean_mod**2d0
                i = i+1
            end do
            corr = 1d0*corr/count
            corr_err = sqrt(1d0*corr_err/count - corr**2d0)
            vel_corr(tau) = corr
            vel_corr_err(tau) = corr_err
            
            mod_corr = 1d0*mod_corr/count
            mod_corr_err = sqrt(1d0*mod_corr_err/count - mod_corr**2d0)
            mod_vel_corr(tau) = mod_corr
            mod_vel_corr_err(tau) = mod_corr_err
            
            !print *, "Tau: ", tau, "# samp: ", count, "Corr: ", corr, "+-", corr_err
        end do

        ! Divide by the correlation at tau=0 to normalize
        ! Propagating the error

        vel_corr_err(:) = abs(vel_corr(:)/vel_corr(0))* &
            sqrt((vel_corr_err(:)/vel_corr(:))**2d0 + (vel_corr_err(0)/vel_corr(0))**2d0)
        vel_corr(:) = vel_corr(:)/vel_corr(0)

        mod_vel_corr_err(:) = mod_vel_corr_err(:)/var_v
        mod_vel_corr(:) = mod_vel_corr(:)/var_v       


    end subroutine calc_veloc_corr

end module