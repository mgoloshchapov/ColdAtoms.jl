using ColdAtoms
using QuantumOptics
using Plots
include("./assets/default.jl")


# # Additional red laser params
red_laser_params = [Ωr, wr, zr];
detuning_params = [Δ0, δ_twophoton(Ωr, Ωb, Δ0)];


# # Period of Rabi oscillations
T0 = T_twophoton(Ωr, Ωb, Δ0);
# # Simulation time
tspan = [0.0:T0/30:2.5*T0;];
# # Initial wavefunction
ψ0 = g;


# # Get initial atom coordinates and velocities
N = 200;
samples, acc_rate = samples_generate(trap_params, atom_params, N; skip=5000, freq=1000);


# # Calculate density matrix and it's squared values to estimate errors
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


# # Take average of state populations over density matrix
Pg = real(expect(g ⊗ dagger(g), ρ_mean)); 
Pp = real(expect(p ⊗ dagger(p), ρ_mean)); 
Pr = real(expect(r ⊗ dagger(r), ρ_mean));


# Plot Rabi oscillations
plot(
    tspan, 
    [Pg Pp Pr], 
    label=["Ground" "Intermediate" "Rydberg"], 
    lw=[3 3 3], 
    color=["blue" "gray" "red"]
    )
title!("State populations")
plot!(legend=:outertop, legendcolumns=3)
xlims!(0.0, maximum(tspan))
ylims!(0.0, 1.05)
xlabel!("Time, μs")
ylabel!("Recapture probability")