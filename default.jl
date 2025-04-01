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

atom_params = [m, kT];
trap_params = [U0, w0, z0];

#Rabi frequencies
Δ0 = 2.0*π * 904.0;
Ωr = 2π * 60.0;
Ωb = 2π * 60.0;
blue_laser_params = [Ωb, wb, zb];
red_laser_params = [Ωr, wr, zr];

detuning_params = [Δ0, δ_twophoton(Ωr, Ωb, Δ0)];
Γ = 2.0*π * 6.0;
decay_params = [Γ/4, 3*Γ/4];