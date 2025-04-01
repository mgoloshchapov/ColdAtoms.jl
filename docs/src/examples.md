# Examples

## Release and recapture experiment

Release and recapture method is used to measure temperature of single neutral atom in optical trap formed by gaussian beam. 

The trap is turned off by a variable amount of time and probability to recapture atom is measured. If we know atom mass and trap parameters, we can extract atom temperature by comparing experimental release and recapture curve with the modelling.

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

# Plot release and recapture curve
figure(figsize=(6,4))
plot(tspan, probs, color="red", linewidth=2)
xlabel(L"Time, $\mu s$")
ylabel("Recapture probability")
gcf()
```