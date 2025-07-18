#Basis states
const basis = NLevelBasis(4);
const g = nlevelstate(basis, 1);
const p = nlevelstate(basis, 2);
const r = nlevelstate(basis, 3);
const gt = nlevelstate(basis, 4);

#Operators
const σgp = g ⊗ dagger(p);
const σpg = p ⊗ dagger(g);
const σpr = p ⊗ dagger(r);
const σrp = r ⊗ dagger(p);
const np = p ⊗ dagger(p);
const nr = r ⊗ dagger(r);
const σgtp = gt ⊗ dagger(p);
const σpgt = p ⊗ dagger(gt);

const σgtr = gt ⊗ dagger(r);
const σrgt = r ⊗ dagger(gt);

const operators = [np, nr, σgp, σpg, σpr, σrp];


#Due to atom dynamics
@inline function Ω(x, y, z, laser_params)
    Ω0, w0, z0, θ, n = laser_params;
    return Ω0 .* A(x, y, z, w0, z0; n=n, θ=θ) .* A_phase(x, y, z, w0, z0; θ=θ);
end;

#Due to Doppler shift for red laser
@inline function Δ(vx, vz, laser_params)
    w0, z0, θ = laser_params[2:4]
    k = 2 * z0/w0^2;

    Δx = k * sin(θ) * vx
    Δz = k * cos(θ) * vz
    return Δx + Δz
end;

#Due to Doppler shifts for red and blue lasers
@inline function δ(vx, vz, red_laser_params, blue_laser_params)
    wr0, zr0, θr = red_laser_params[2:4];
    wb0, zb0, θb = blue_laser_params[2:4];
    
    kr = 2 * zr0/wr0^2;
    kb = 2 * zb0/wb0^2;

    δx = (kr*sin(θr) + kb*sin(θb))*vx
    δz = (kr*cos(θr) + kb*cos(θb))*vz 

    return δx + δz
end;


#Two-photon Rydberg hamiltonian for 1 atom
function Hamiltonian(Ωr, Ωb, Δ, δ)
    return TimeDependentSum(
        [
            t -> -Δ(t),
            t -> -δ(t),
            t -> Ωr(t) ./2.0,
            t -> conj.(Ωr(t)) ./2.0,
            t -> Ωb(t)/2.0,
            t -> conj.(Ωb(t)) ./2.0,
        ],
        
        [
            np,
            nr,
            σgp,
            σpg,
            σpr,
            σrp  
        ]
    )
end;



#Jump operators for master equation 
@inline function JumpOperators(decay_params)
    Γg, Γgt, Γr = decay_params;
    return [[sqrt(Γg)*σgp, sqrt(Γgt)*σgtp, sqrt(Γr)*σgtr], 
    [sqrt(Γg)*σpg, sqrt(Γgt)*σpgt, sqrt(Γr)*σrgt]]
end;


@inline function GenerateHamiltonian(
    sample, 
    ωr, ωz,
    free_motion,
    atom_motion,

    tspan_noise,
    f,
    red_laser_phase_amplitudes,
    blue_laser_phase_amplitudes,
    nodes,

    red_laser_params,
    blue_laser_params,

    Δ0, 
    δ0
    )
    # Trajectories
    xi, yi, zi, vxi, vyi, vzi = atom_motion ? sample : zeros(6);
    X  = t -> R(t, xi, vxi, ωr; free=free_motion);
    Y  = t -> R(t, yi, vyi, ωr; free=free_motion);
    Z  = t -> R(t, zi, vzi, ωz; free=free_motion);
    Vx = t -> V(t, xi, vxi, ωr; free=free_motion);
    Vz = t -> V(t, zi, vzi, ωz; free=free_motion);

    # Generate phase noise traces for red and blue lasers
    ϕ_red_res  = ϕ(tspan_noise, f, red_laser_phase_amplitudes);
    ϕ_blue_res = ϕ(tspan_noise, f, blue_laser_phase_amplitudes);

    # Interpolate phase noise traces to pass to hamiltonian
    ϕ_red  = interpolate(nodes, ϕ_red_res, Gridded(Linear()));
    ϕ_blue = interpolate(nodes, ϕ_blue_res, Gridded(Linear()));


    # Hamiltonian params trajectories
    Ωr = t -> exp(1.0im * ϕ_red(t)) * Ω(X(t), Y(t), Z(t), red_laser_params);
    Ωb = t -> exp(1.0im * ϕ_blue(t)) * Ω(X(t), Y(t), Z(t), blue_laser_params);

    Ht = TimeDependentSum(
        [
            t -> -Δ(Vx(t), Vz(t), red_laser_params) - Δ0,
            t -> -δ(Vx(t), Vz(t), red_laser_params, blue_laser_params) - δ0,
            t -> Ωr(t) / 2.0,
            t -> conj(Ωr(t)) / 2.0,
            t -> Ωb(t) / 2.0,
            t -> conj(Ωb(t)) / 2.0,
        ],
        operators
        );

    return Ht
end;


mutable struct RydbergConfig
    tspan::Vector{Float64}
    ψ0::Ket{NLevelBasis{Int64}, Vector{ComplexF64}}

    atom_params::Vector{Float64}
    trap_params::Vector{Float64}
    n_samples::Int64

    f::Vector{Float64}
    red_laser_phase_amplitudes::Vector{Float64}
    blue_laser_phase_amplitudes::Vector{Float64}
    
    red_laser_params::Vector{Float64}
    blue_laser_params::Vector{Float64}
    
    detuning_params::Vector{Float64}
    decay_params::Vector{Float64}
    atom_motion::Bool
    free_motion::Bool
    laser_noise::Bool
    spontaneous_decay_intermediate::Bool
    spontaneous_decay_rydberg::Bool
end


function simulation(cfg::RydbergConfig)
    samples = samples_generate(
        cfg.trap_params,
        cfg.atom_params,
        cfg.n_samples;
        harmonic=false
        )[1]

    ωr, ωz = trap_frequencies(cfg.atom_params, cfg.trap_params);
    Δ0, δ0 = cfg.detuning_params;

    Γg, Γgt = cfg.spontaneous_decay_intermediate ? cfg.decay_params[1:2] : [0.0, 0.0]
    Γr      = cfg.spontaneous_decay_rydberg      ? cfg.decay_params[3]   :  0.0
    decay_params = [Γg, Γgt, Γr]
    J, Jdagger = JumpOperators(decay_params)

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

    for sample in ProgressBars.ProgressBar(samples)
        Ht = GenerateHamiltonian(
            sample, 
            ωr, ωz,
            cfg.free_motion,
            cfg.atom_motion,
        
            tspan_noise,
            cfg.f,
            red_laser_phase_amplitudes,
            blue_laser_phase_amplitudes,
            nodes,
        
            cfg.red_laser_params,
            cfg.blue_laser_params,
        
            Δ0, 
            δ0
            )

        super_operator(t, rho) = Ht, J, Jdagger
        _, ρ_temp = timeevolution.master_dynamic(cfg.tspan, ρ0, super_operator);

        ρ_mean  .+= ρ_temp
        ρ2_mean .+= ρ_temp .^ 2
    end;

    return ρ_mean ./ cfg.n_samples, ρ2_mean ./ cfg.n_samples
end;

function Ω_twophoton(Ωr, Ωb, Δ)
    return Ωb * Ωr / (2.0 * Δ)
end;

function T_twophoton(Ωr, Ωb, Δ)
    return 2.0*π / Ω_twophoton(Ωr, Ωb, Δ)
end;

function δ_twophoton(Ωr, Ωb, Δ)
    return (Ωb^2 - Ωr^2)/(4.0 * Δ)
end;

function Ωr_required(Ω, Ωb, Δ)
    return 2.0 * Δ * Ω / Ωb
end;
