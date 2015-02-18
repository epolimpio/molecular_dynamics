program ArgonSim

    use simulation
    use simParam
    
    call init_all
    call stabilize_temp
    call simulate_dynamics

end program ArgonSim