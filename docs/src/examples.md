# Examples

## Release and recapture experiment

```@example
using ColdAtoms
using PyPlot

# Atom and trap parameters
U0, w0, λ, M2 = 1000.0, 1.1, 0.852, 1.3;
z0 = ColdAtoms.w0_to_z0(w0, λ, M2);
m, T = 86.9091835, 50.0;
trap_params = [U0, w0, z0];
atom_params = [m, T];

# Simulation parameters
tspan = [0.0:1.0:40.0;];
N = 10000 

# Calculate recapture probabilities
probs, acc_rate = ColdAtoms.release_recapture(
    tspan, 
    trap_params, 
    atom_params, 
    N; 
    harmonic=false);

figure()
plot(tspan, probs)
```