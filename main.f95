program ArgonSim

    use verlet

    call init_all
    call stabilize_temp
    call simulate_dynamics

end program ArgonSim