using Distributed

function parallel_evolution(smpl,cnfg,
    tspan_noise, red_laser_phase_amplitudes,
    blue_laser_phase_amplitudes, nodes)
    
    ωr, ωz = trap_frequencies(cnfg.atom_params, cfg.trap_params);
    Δ0, δ0 = cnfg.detuning_params;

    Γ0, Γ1, Γl   = cfg.spontaneous_decay_intermediate ? cfg.decay_params[1:3] : zeros(3)
    Γr           = cfg.spontaneous_decay_rydberg      ? cfg.decay_params[4]   :  0.0
    decay_params = [Γ0, Γ1, Γl, Γr]
    J, Jdagger = JumpOperators(decay_params)

    ρ_0 = cfg.ψ0 ⊗ dagger(cfg.ψ0);

    Ht = GenerateHamiltonian(
        smpl,  ωr, ωz,
        cnfg.free_motion, cnfg.atom_motion,
        tspan_noise, cnfg.f,
        red_laser_phase_amplitudes,
        blue_laser_phase_amplitudes,
        nodes,         
        cnfg.red_laser_params,
        cnfg.blue_laser_params,
        Δ0, δ0        )

    super_operator(t, rho) = Ht, J, Jdagger
    _, ρ_temp = timeevolution.master_dynamic(cnfg.tspan, ρ_0, super_operator);
    return ρ_temp
end


"""    ωr, ωz = trap_frequencies(cfg.atom_params, cfg.trap_params);
    Δ0, δ0 = cfg.detuning_params;
    Γ0, Γ1, Γl   = cfg.spontaneous_decay_intermediate ? cfg.decay_params[1:3] : zeros(3)
    Γr           = cfg.spontaneous_decay_rydberg      ? cfg.decay_params[4]   :  0.0
    decay_params = [Γ0, Γ1, Γl, Γr]
    J, Jdagger = JumpOperators(decay_params)"""

function simulation_parallel(cfg::RydbergConfig)
    samples = samples_generate(
        cfg.trap_params,
        cfg.atom_params,
        cfg.n_samples;
        harmonic=false
        )[1]

    #ωr, ωz = trap_frequencies(cfg.atom_params, cfg.trap_params);
    #Δ0, δ0 = cfg.detuning_params;
   """ Γ0, Γ1, Γl   = cfg.spontaneous_decay_intermediate ? cfg.decay_params[1:3] : zeros(3)
    Γr           = cfg.spontaneous_decay_rydberg      ? cfg.decay_params[4]   :  0.0
    decay_params = [Γ0, Γ1, Γl, Γr]
    J, Jdagger = JumpOperators(decay_params)"""

    ρ0 = cfg.ψ0 ⊗ dagger(cfg.ψ0);
    #Density matrix averaged over realizations of laser noise and atom dynamics.
    ρ_mean  = [zero(ρ0) for _ ∈ 1:length(cfg.tspan)];
    ρ_temp  = [zero(ρ0) for _ ∈ 1:length(cfg.tspan)];
    #Second moment for error estimation of level populations. 
    ρ2_mean = [zero(ρ0) for _ ∈ 1:length(cfg.tspan)];

    tspan_noise = [0.0:cfg.tspan[end]/1000:cfg.tspan[end];];
    nodes = (tspan_noise, );
    red_laser_phase_amplitudes  = cfg.laser_noise ? cfg.red_laser_phase_amplitudes  : zero(cfg.red_laser_phase_amplitudes);
    blue_laser_phase_amplitudes = cfg.laser_noise ? cfg.blue_laser_phase_amplitudes : zero(cfg.blue_laser_phase_amplitudes);

    #addprocs(4) #num_procs)  # Adjust the number of workers as needed

    results = pmap(smpl -> parallel_evolution(smpl,cnfg,tspan_noise,red_laser_phase_amplitudes,blue_laser_phase_amplitudes, nodes), samples)
    
    for res in results
        ρ_mean  .+= res
        ρ2_mean .+= res .^ 2
    end;

    rmprocs(workers())

    return ρ_mean ./ cfg.n_samples, ρ2_mean ./ cfg.n_samples
end;
