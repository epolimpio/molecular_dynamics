program ArgonSim

    use simulation
    use output_file
    use simParam
    
    real(8), allocatable :: settings(:,:)
    integer :: i
    character(len=3) :: str_file

    call init_temps
    do i = 1, size(settings, 2)
        call set_param(864, settings(1,i), settings(2,i))
        write(str_file, "(I0.3)") i
        call init_all(.FALSE.)
        call init_files(str_file)
        call stabilize_temp
        call simulate_dynamics
        call close_files
    end do

contains
    
    subroutine init_temps

        integer :: i, j, count
        real(8), parameter :: temps(1:8) = (/ 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0 /)
        real(8), parameter :: dens(1:8) = (/ 0.3, 0.45, 0.6, 0.75, 0.88, 0.95, 1.0, 1.2 /)

        allocate(settings(2,size(temps)*size(dens)))
        count = 0
        do i = 1, size(temps)
            do j = 1, size(dens)
                count = count + 1
                settings(1,count) = temps(i)
                settings(2,count) = dens(j)
            end do
        end do

    end subroutine init_temps

end program ArgonSim