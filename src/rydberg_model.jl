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

const operators = [np, nr, σgp, σpg, σpr, σrp];


#Due to atom dynamics
function Ω(x, y, z, laser_params; n=1, θ=0.0)
    Ω0, w0, z0 = laser_params;
    return Ω0 .* A(x, y, z, w0, z0; n=n, θ=θ) .* A_phase(x, y, z, w0, z0; θ=θ);
end;

#Due to Doppler shift for red laser
function Δ(vx, vz, laser_params; θ=0.0)
    Ω0, w0, z0 = laser_params
    k = 2 * z0/w0^2;

    Δx = k * sin(θ) * vx
    Δz = k * cos(θ) * vz
    return Δx + Δz
end;

#Due to Doppler shifts for red and blue lasers
function δ(vx, vz, red_laser_params, blue_laser_params; θr=0.0, θb=0.0)
    Ωr0, wr0, zr0 = red_laser_params;
    Ωb0, wb0, zb0 = blue_laser_params;
    
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
function JumpOperators(decay_params)
    Γg, Γgt = decay_params;
    return [sqrt(Γg)*σgp, sqrt(Γgt)*σgtp], [sqrt(Γg)*σpg, sqrt(Γgt)*σpgt]
end;



"""
    simulation(
        tspan, ψ0, 
        atom_params,
        trap_params,
        samples,
        f,
        red_laser_phase_amplitudes,
        blue_laser_phase_amplitudes,
        red_laser_params,
        blue_laser_params,
        detuning_params,
        decay_params;
        atom_motion=true,
        free_motion=true,
        laser_noise=true,
        spontaneous_decay=true,
        parallel=false
        )

Simulate two-photon Rydberg excitation of single atom with several sources of decoherence

### Input

- `tspan` -- vector specifying the points of time for which output should be displayed
- `ψ0` -- initial wavefunction vector of normalized complex amplitudes ``[c_{g}, c_{p}, c_{r}]``
- `atom_params` -- vector [atom mass in a.u., atom temperature in ``\\mu K``]
- `trap_params` -- vector [trap depth ``U_{0}`` in ``\\mu K``, beam waist radius in ``\\mu m``, beam Rayleigh length in ``\\mu m``]
- `samples` -- Monte-Carlo samples of initial atom coordinates and velocities, can be received using samples_generate
- `f` -- frequencies at which laser phase noise is sampled
- `red_laser_phase_amplitudes` -- amplitudes of red laser phase noise for correspoding frequencies `f`
- `blue_laser_phase_amplitudes` -- amplitudes of blue laser phase noise for correspoding frequencies `f`
- `red_laser_params` -- write explanation
- `blue_laser_params` -- write explanation
- `detuning_params` -- vector [Δ0, δ0], which sets detuning from intermediate level and Rydberg level correspondingly
- `decay_params` -- write explanation
- `atom_motion` -- (optional, default: `true`) if set to true, atom motion is included
- `free_motion` -- (optional, default: `true`) if set to true, trap is turned off
- `laser_noise` -- (optional, default: `true`) if set to true, laser phase noise is included
- `spontaneous_decay` -- (optional, default: `true`) if set to true, spontaneous decay from intermediate level is included
- `parallel` -- (optional, default: `false`) parallel implementation is under development
- `n` -- (optional, default: `1`) super-gauss parameter


### Output

Outputs Monte-Carlo averaged density matrix and squared density matrix for error calculation
with elements ordered in correspondence with order ground, intermediate, Rydberg

"""
function simulation(
    tspan, ψ0, 
        
    atom_params,
    trap_params,
    samples,
    
    f,
    red_laser_phase_amplitudes,
    blue_laser_phase_amplitudes,
    
    red_laser_params,
    blue_laser_params,
    
    detuning_params,
    decay_params;

    atom_motion=true,
    free_motion=true,
    laser_noise=true,
    spontaneous_decay=true,
    θr=0.0,
    θb=0.0,
    n=1
    )
    N = length(samples);

    ωr, ωz = trap_frequencies(atom_params, trap_params);
    Δ0, δ0 = detuning_params;

    if spontaneous_decay
        decay_params_temp = decay_params;
    else
        decay_params_temp = [0.0, 0.0];
    end;

    Γg, Γgt = decay_params_temp;
    J, Jdagger = [sqrt(Γg)*σgp, sqrt(Γgt)*σgtp], [sqrt(Γg)*σpg, sqrt(Γgt)*σpgt];

    ρ0 = ψ0 ⊗ dagger(ψ0);

    #Density matrix averaged over realizations of laser noise and atom dynamics.
    ρ_mean = [zero(ψ0 ⊗ dagger(ψ0)) for _ ∈ 1:length(tspan)];
    ρ_temp = [zero(ψ0 ⊗ dagger(ψ0)) for _ ∈ 1:length(tspan)];

    #Second moment for error estimation of level populations. 
    ρ2_mean = [zero(ψ0 ⊗ dagger(ψ0)) for _ ∈ 1:length(tspan)];

    tspan_noise = [0.0:tspan[end]/1000:tspan[end];];
    nodes = (tspan_noise, );

    if laser_noise
        red_laser_phase_amplitudes_temp = red_laser_phase_amplitudes;
        blue_laser_phase_amplitudes_temp = blue_laser_phase_amplitudes;
    else
        red_laser_phase_amplitudes_temp = zero(red_laser_phase_amplitudes);
        blue_laser_phase_amplitudes_temp = zero(blue_laser_phase_amplitudes);
    end;

    
    for i ∈ 1:N
        if atom_motion
            #Atom initial conditions
            xi, yi, zi, vxi, vyi, vzi = samples[i];
        else
            xi, yi, zi, vxi, vyi, vzi = zeros(6);
        end;
        
        #Atom trajectories
        X = t -> R(t, xi, vxi, ωr; free=free_motion);
        Y = t -> R(t, yi, vyi, ωr; free=free_motion);
        Z = t -> R(t, zi, vzi, ωz; free=free_motion);
        Vx = t -> V(t, xi, vxi, ωr; free=free_motion);
        Vz = t -> V(t, zi, vzi, ωz; free=free_motion);
        
 
        #Generate phase noise traces for red and blue lasers
        ϕ_red_res = ϕ(tspan_noise, f, red_laser_phase_amplitudes_temp);
        ϕ_blue_res = ϕ(tspan_noise, f, blue_laser_phase_amplitudes_temp);

        #Interpolate phase noise traces to pass to hamiltonian
        ϕ_red = interpolate(nodes, ϕ_red_res, Gridded(Linear()));
        ϕ_blue = interpolate(nodes, ϕ_blue_res, Gridded(Linear()));

        #Hamiltonian params trajectories
        Ht = TimeDependentSum(
        [
            t -> -Δ(Vx(t), Vz(t), red_laser_params; θ=θr) - Δ0,
            t -> -δ(Vx(t), Vz(t), red_laser_params, blue_laser_params; θr=θr, θb=θb) - δ0,
            t -> exp(1.0im * ϕ_red(t)) * Ω(X(t), Y(t), Z(t), red_laser_params; n=n, θ=θr) / 2.0,
            t -> conj(exp(1.0im * ϕ_red(t)) * Ω(X(t), Y(t), Z(t), red_laser_params; n=n, θ=θr) / 2.0),
            t -> exp(1.0im * ϕ_blue(t)) * Ω(X(t), Y(t), Z(t), blue_laser_params; n=n, θ=θb) / 2.0,
            t -> conj(exp(1.0im * ϕ_blue(t)) * Ω(X(t), Y(t), Z(t), blue_laser_params; n=n, θ=θb) / 2.0)
        ],
        operators
        );

        # #Returns hamiltonian and jump operators in a form required by timeevolution.master_dynamic
        function super_operator(t, rho)
            return Ht, J, Jdagger;
        end;
        
        _, ρ_temp = timeevolution.master_dynamic(tspan, ρ0, super_operator);

        ρ_mean = ρ_mean + ρ_temp;
        ρ2_mean = ρ2_mean + ρ_temp .^ 2;
    end;

    return ρ_mean/N, ρ2_mean/N
end;

function Ω_twophoton(Ωr, Ωb, Δ)
    return Ωb * Ωr / (2.0 * Δ)
end;

function T_twophoton(Ωr, Ωb, Δ)
    return 2.0*π / Ω_twophoton(Ωr, Ωb, Δ);
end;

function δ_twophoton(Ωr, Ωb, Δ)
    return (Ωb^2 - Ωr^2)/(4.0 * Δ)
end;

function Ωr_required(Ω, Ωb, Δ)
    return 2.0 * Δ * Ω / Ωb
end;
