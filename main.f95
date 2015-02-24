program ArgonSim

    use simulation
    use output_file
    use simParam
    
    call init_all(.TRUE.)
    call init_files("Liquid")
    call stabilize_temp
    call simulate_dynamics
    call close_files

end program ArgonSim