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


#Rb87 mass in a.u.
m = 86.9091835;       

#Params for laser phase noise
h0 = 13.0 * 1e-6;    #MHz^2/MHz
hg1 = 25.0 * 1e-6;   #MHz^2/MHz
hg2 = 10.0e3 * 1e-6; #MHz^2/MHz
fg1 = 130.0 * 1e-3;  #MHz
fg2 = 234.0 * 1e-3;  #MHz
σg1 = 18.0 * 1e-3;   #MHz
σg2 = 1.5 * 1e-3;    #MHz

red_laser_phase_params  = [h0, [hg1, hg2], [σg1, σg2], [fg1, fg2]];
blue_laser_phase_params = [h0, [hg1, hg2], [σg1, σg2], [fg1, fg2]];
f = [0.01:0.0025:1.0;];
red_laser_phase_amplitudes = ϕ_amplitudes(f, red_laser_phase_params);
blue_laser_phase_amplitudes = ϕ_amplitudes(f, blue_laser_phase_params);

#Excitation beam parameters
λr = 0.795;
λb = 0.475;
wr = 10.0;
wb = 3.5;
zr = w0_to_z0(wr, λr);
zb = w0_to_z0(wr, λb);

# Atom and trap parameters
T = 50.0;
U0 = 1000.0;
w0 = 1.1;
λ0 = 0.813;
M2 = 1.3;
z0 = w0_to_z0(w0, λ0, M2);
atom_params = [m, T];
trap_params = [U0, w0, z0];

#Rabi frequencies
Δ0 = 2.0*π * 904.0;
Ωr = 2π * 60.0;
Ωb = 2π * 60.0;
blue_laser_params = [Ωb, wb, zb];
red_laser_params = [Ωr, wr, zr];

# Detunings and decay params
detuning_params = [Δ0, δ_twophoton(Ωr, Ωb, Δ0)];
Γ = 2.0*π * 6.0;
decay_params = [Γ/4, 3*Γ/4];

# Period of Rabi oscillations
T0 = T_twophoton(Ωr, Ωb, Δ0);
# Simulation time
tspan = [0.0:T0/30:2.5*T0;];
# Initial wavefunction
ψ0 = g;



# # Get initial atom coordinates and velocities
N = 200;
samples, acc_rate = samples_generate(trap_params, atom_params, N; skip=5000, freq=1000);



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
plot(
    tspan, 
    [Pg Pp Pr], 
    label=["Ground" "Intermediate" "Rydberg"], 
    width=[3 3 3], 
    color=["blue" "gray" "red"]
    )
plot!(legend=:outertop, legendcolumns=3)
title!("State populations")
xlims!(0.0, maximum(tspan))
ylims!(0.0, 1.05)
xlabel!("Time, μs")
ylabel!("Probability")
```