using ColdAtoms
using QuantumOptics


m = 86.9091835;
U0 = 340.0;
T = 30.0;
w0 = 1.2;
λ0 = 0.810;
M2 = 1.0;
z0 = w0_to_z0(w0, λ0, M2);
atom_params = [m, T];
trap_params = [U0, w0, z0];


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
λb = 0.474;
wr = 40.0;
wb = 5.0;
zr = w0_to_z0(wr, λr);
zb = w0_to_z0(wr, λb);

#Rabi frequencies
Δ0 = 2.0*π * 904.0;
Ωr = 2π * 80.0;
Ωb = 2π * 80.0;
blue_laser_params = [Ωb, wb, zb];
red_laser_params = [Ωr, wr, zr];


# Detunings and decay params
detuning_params = [Δ0, δ_twophoton(Ωr, Ωb, Δ0)];
Γ = 2.0*π * 6.0;
decay_params = [Γ/4, 3*Γ/4];


# Period of Rabi oscillations
T0 = T_twophoton(Ωr, Ωb, Δ0);
tspan = [0.0:T0/20:3*T0;]; 

# Get initial atom coordinates and velocities
ψ0 = g;


N = 100;
samples, acc_rate = samples_generate(trap_params, [m, t] , N; skip=5000, freq=1000);

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

        spontaneous_decay=false,
        atom_motion=true,
        laser_noise=false,
        θr=0.0,
        θb=π/2,
        free_motion=false,
    );

# Take average of state populations over density matrix
Pg = real(expect(g ⊗ dagger(g), ρ_mean));
Pp = real(expect(p ⊗ dagger(p), ρ_mean));
Pr = real(expect(r ⊗ dagger(r), ρ_mean));
