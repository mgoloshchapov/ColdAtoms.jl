# Examples

## Release and recapture experiment

Release and recapture method is used to measure temperature of single neutral atom in optical trap formed by gaussian beam. 

The trap is turned off by a variable amount of time and probability to recapture atom is measured. If we know atom mass and trap parameters, we can extract atom temperature by comparing experimental release and recapture curve with the modelling.

```@example
using ColdAtoms
using Plots

# Atom and trap parameters
U0, w0, λ, M2 = 1000.0, 1.1, 0.852, 1.3;
z0 = w0_to_z0(w0, λ, M2);
m, T = 86.9091835, 50.0;
trap_params = [U0, w0, z0];
atom_params = [m, T];

# Simulation parameters
tspan = [0.0:1.0:50.0;];
N = 10000 

# Calculate recapture probabilities
probs, acc_rate = release_recapture(
    tspan, 
    trap_params, 
    atom_params, 
    N; 
    harmonic=false);

# Plot release and recapture curve
plot(tspan, probs, label=false, width=3, color="red")
xlims!(0.0, 50.0)
ylims!(0.0, 1.05)
xlabel!("Time, μs")
ylabel!("Recapture probability")
```

## Two-photon Rydberg excitation 

Two-photon Rydberg excitation is used, for example, to implement native CZ gate on neutral atom quantum computers. Here we implement modelling of two-photon Rydberg excitation with different sources of decoherence: atom dynamics, laser phase noise, spontaneous decay from intermediate state. 

```@example
using ColdAtoms
using QuantumOptics
using Plots
include("./assets/default.jl")


# Additional red laser params
red_laser_params = [Ωr, wr, zr];
detuning_params = [Δ0, δ_twophoton(Ωr, Ωb, Δ0)];


# Period of Rabi oscillations
T0 = T_twophoton(Ωr, Ωb, Δ0);
# Simulation time
tspan = [0.0:T0/30:2.5*T0;];
# Initial wavefunction
ψ0 = g;


# Calculate density matrix and it's squared values to estimate errors
ρ_mean, ρ2_mean = 
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
        
        spontaneous_decay=true,
        atom_motion=true,
        laser_noise=true,
        parallel=false,
        free_motion=true
    );


# Take average of state populations over density matrix
Pg = real(expect(g ⊗ dagger(g), ρ_mean)); 
Pp = real(expect(p ⊗ dagger(p), ρ_mean)); 
Pr = real(expect(r ⊗ dagger(r), ρ_mean));


# Plot Rabi oscillations
plot(tspan, [Pg, Pp, Pr], label=["Ground population", "Intermediate population", "Rydberg population"], width=3, color=["blue", "orange", "red"])
plot!(legend=:outerbottom, legendcolumns=3)
xlims!(0.0, maximum(tspan))
ylims!(0.0, 1.05)
xlabel!("Time, μs")
ylabel!("Recapture probability")
```