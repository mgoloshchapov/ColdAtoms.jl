function is_zero(x)
    return x == 0
end;


function release_evolve(tspan, cord, atom_params, trap_params; eps=1e-3)
    xi, yi, zi, vxi, vyi, vzi = cord;
    m, T = atom_params;
    U0, w0, z0 = trap_params;
    
    x = xi .+ vxi * tspan;
    y = yi .+ vyi * tspan - g0 * tspan .^2;
    z = zi .+ vzi * tspan;
    
    kinetic = K(cord, trap_params, m);
    potential = U0 .* (1.0 .- A(x, y, z, w0, z0) .^2);
    recap = (kinetic .+ potential) .< U0 * (1.0-eps);
    
    idx = findfirst(is_zero, recap);
    
    #Changed != nothing to !isnothing
    if !isnothing(idx)
        recap[idx:end] .= 0;
    end;
    
    return recap
end; 


"""
    release_recapture(tspan, trap_params, atom_params, N; freq=10, skip=1000, eps=1e-3, harmonic=true)

Simulate release and recapture experiment to estimate atom's temperature.

### Input

- `tspan` -- vector specifying the points of time for which output should be displayed
- `trap_params` -- vector [trap depth U0 in ``\\mu K``, beam waist radius in ``\\mu m``, beam Rayleigh length in ``\\mu m``]
- `atom_params` -- vector [atom mass in a.u., atom temperature in ``\\mu K``]
- `N` -- number of Monte-Carlo samples, the same as number of atoms
- `freq` -- (optional, default: `10`) number of Metropolis steps skipped between samples to reduce sample dependency
- `skip` -- (optional, default: `1000`) number of Metropolis steps skipped before the Markov Chain is considered to reach stationary distribution
- `eps` -- (optional, default: `1e-3`) cutoff to regularize Metropolis sampler, atoms that have energy over ``U_{0}(1-eps)`` are considered to be out of trap
- `harmonic` -- (optional, default: `true`) uses harmonic approximation of gaussian beam if set to `true`, otherwise uses Metropolis sampler

### Output

List of recapture probabilities corresponding to times in `tspan` and acceptance rate of Metropolis algorithm. 
If `harmonic` is set to `true`, acceptance rate is set to 1.0. 
"""
function release_recapture(tspan, trap_params, atom_params, N; freq=10, skip=1000, eps=1e-3, harmonic=true)
    samples, acc_rate = samples_generate(trap_params, atom_params, N; freq=freq, skip=skip, harmonic=harmonic);
    recapture = zeros(length(tspan));
    
    for i ∈ 1:N
        recapture += release_evolve(tspan, samples[i], atom_params, trap_params; eps=eps);
    end;
    
    return recapture ./ N, acc_rate
end;