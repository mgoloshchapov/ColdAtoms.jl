basis_fidelity_states = [
    ket_0, 
    ket_1,
    (ket_0 + ket_1)/sqrt(2),
    (ket_0 - ket_1)/sqrt(2),
    (ket_0 + 1.0im * ket_1)/sqrt(2),
    (ket_0 - 1.0im * ket_1)/sqrt(2)
    ]


function get_rydberg_fidelity_configs(cfg, n_samples=20)
    configs = OrderedDict()

    # Config to measure error from intermediate state decay
    cfg_t = deepcopy(cfg)
    cfg_t.atom_params[2] = 1.0
    cfg_t.spontaneous_decay_intermediate = true
    cfg_t.spontaneous_decay_rydberg      = false
    cfg_t.laser_noise = false
    cfg_t.free_motion = true
    cfg_t.n_samples = 1
    configs["Intermdeiate state decay"] = cfg_t

    # Config to measure error from rydberg state decay
    cfg_t = deepcopy(cfg)
    cfg_t.atom_params[2] = 1.0
    cfg_t.spontaneous_decay_intermediate = false
    cfg_t.spontaneous_decay_rydberg      = true
    cfg_t.laser_noise = false
    cfg_t.free_motion = true
    cfg_t.n_samples = 1
    configs["Rydberg state decay"] = cfg_t

    # Config to measure error from laser_noise
    cfg_t = deepcopy(cfg)
    cfg_t.atom_params[2] = 1.0
    cfg_t.spontaneous_decay_intermediate = false
    cfg_t.spontaneous_decay_rydberg = false
    cfg_t.laser_noise = true
    cfg_t.free_motion = false
    cfg_t.n_samples = n_samples
    configs["Laser noise"] = cfg_t

    # Config to measure error from temperature
    cfg_t = deepcopy(cfg)
    cfg_t.spontaneous_decay_intermediate = false
    cfg_t.spontaneous_decay_rydberg = false
    cfg_t.laser_noise = false
    cfg_t.free_motion = true
    cfg_t.n_samples = n_samples
    configs["Atom motion"] = cfg_t


    # Config to measure total error
    cfg_t = deepcopy(cfg)
    cfg_t.spontaneous_decay_intermediate = true
    cfg_t.spontaneous_decay_rydberg = true
    cfg_t.laser_noise = true
    cfg_t.free_motion = true
    cfg_t.n_samples = n_samples
    configs["Total"] = cfg_t

    return configs
end 


function get_rydberg_infidelity(
    cfg::RydbergConfig;
    U=dense(identityoperator(ColdAtoms.basis)), 
    states=basis_fidelity_states, 
    n_samples=1,
    ode_kwargs...)

    configs = get_rydberg_fidelity_configs(cfg)
    names = collect(keys(configs))
    infidelities = Dict()

    for name in ProgressBar(names)
        cfg_t = deepcopy(configs[name])
        println("Measuring error from $(name)...")
        infidelity_avg = 0.0
        for state in states
            ψ_ideal = U * state;
            cfg_t.ψ0 = state
            ρ_real = simulation(cfg_t)[1][end]
            infidelity_avg += 1.0 - real(dagger(ψ_ideal) * ρ_real * ψ_ideal)
        end
        infidelities[name] = infidelity_avg / length(states)

        println()
        println("Infidelity from $(name): $(round(100.0*infidelities[name]; digits=4)) %")
    end

    return infidelities
end

function plot_rydberg_infidelity(infidelities)
    keys_iF   = collect(keys(infidelities))
    values_iF = 100 .* collect(values(infidelities))

    bar(keys_iF, values_iF; 
    ylabel="Infidelity, %", 
    title="Error budget for 10π pulse",
    angle=45,
    label=nothing)
end

function get_parity_fidelity(cfg::CZLPConfig; ode_kwargs...)
    cfg_parity = deepcopy(cfg)
    ket_pos = (ket_0 + ket_1) / sqrt(2)
    cfg_parity.ψ0 = ket_pos ⊗ ket_pos

    ρ = simulation_czlp(cfg_parity; ode_kwargs...)[1][end]

    Had = Id ⊗ ColdAtoms.Hadamard
    ρ .=  Had * ρ * dagger(Had)
    Phi_p = (ket_0 ⊗ ket_0 + ket_1 ⊗ ket_1)/sqrt(2)

    F = real(dagger(Phi_p) * ρ * Phi_p)
    return F, ρ
end


function get_parity_fidelity_temp(ρ, ϕ1; ode_kwargs...)
    Had = Id ⊗ ColdAtoms.Hadamard
    global_RZ = ColdAtoms.RZ(ϕ1) ⊗ ColdAtoms.RZ(ϕ1);

    ρt =  Had * global_RZ * ρ * dagger(global_RZ) * dagger(Had)
    Phi_p = (ket_0 ⊗ ket_0 + ket_1 ⊗ ket_1)/sqrt(2)

    F = real(dagger(Phi_p) * ρt * Phi_p)
    return F, ρt
end