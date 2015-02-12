program ArgonSim

    use initial_conditions
    use verlet
    use physical_quantities
    use plot_part
    use simParam
    !use sub_regions
    implicit none

    ! Initialize array structure
    real(8) :: Ek, Ep, temp_out
    real(8) :: mom(3)
    integer(8) :: i, frame
    character(5) :: strframe
    logical :: renorm
    integer :: renorm_count
    
    ! Counter frame used to generate positions data
    frame = 0

    ! Start with renormalization
    renorm = .TRUE.
    renorm_count = 0
    


    ! Initialize random distributions & lattices
    
    
    ! Initialize graphics
    !call InitPlot('lightblue', 700, 700, 'out', 1)
    !call Framing(0._8, 0._8, l, l)
    !call PutStopButton()

    ! Calculate initial energies and initialize accelerations
    call kinetic_energy(v, Ek, mom, temp_out)
    call total_force(l, r, a, Ep)

    open(unit = 1001, file = "data/energy.dat")
    open(unit = 1002, file = "data/temperature.dat")

    ! Print initial data
    print *, "Initial Momenta: ", mom
    print *, "Initial Ek: ", Ek
    print *, "Initial Ep: ", Ep
    print *, "Initial E: ", Ek+Ep
    do i = 1, NUM_STEPS
        ! Runs the verlet time evolution for each time step
        call verlet_algorithm(l, DT, r, v, a, Ep)
        call kinetic_energy(v, Ek, mom, temp_out)
        write(1001, *) Ek, Ep, Ek+Ep, temp_out
        call plot_particles(r, N)
        ! Check for renormalization and renormalize
        if (mod(i, STEPS_TO_RENORM) == 0 .AND. renorm) then
            print *, "Renorm"
            renorm_count = renorm_count + 1
            v = sqrt(1d0*TEMP/temp_out)*v
            ! Check if temperature has drift much
            if (renorm_count == NUM_RENORM) then
                renorm = .FALSE.
            end if
        end if
        ! Monitor the results printing positions and energies to files
        if (mod(i,100) == 0) then
            frame = frame + 1
            print *, frame, "Ek: ", Ek, "Ep: ", Ep, "E: ", Ek+Ep, "T: ", temp_out
            write(strframe,'(I5)') frame
        end if

    end do

    !call EndPlot()
    close(1001)

end program ArgonSim

