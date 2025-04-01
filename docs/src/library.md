# ColdAtoms.jl Documentation

```@docs
release_recapture(tspan, trap_params, atom_params, N; <keyword arguments>)
```

```@docs
w0_to_z0(w0, λ, M2)
```

```@docs
trap_frequencies(atom_params, trap_params)
```

```@docs
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

        <keyword arguments>
        )
```
