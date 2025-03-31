"""
trap_params: [U0(μK), w0(μm), z0(μm)]
atom_params: [m(a.u.), T(μK)], T is used as initial guess for temperature
tspan      : time from experiment(in μs)
probability: probability to recapture atom from experiment
"""

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


function release_recapture(tspan, trap_params, atom_params, N; freq=10, skip=1000, eps=1e-3, harmonic=true)
    samples, acc_rate = samples_generate(trap_params, atom_params, N; freq=freq, skip=skip, harmonic=harmonic);
    recapture = zeros(length(tspan));
    
    for i ∈ 1:N
        recapture += release_evolve(tspan, samples[i], atom_params, trap_params; eps=eps);
    end;
    
    return recapture ./ N, acc_rate
end;
# function release_recapture_antitrapping(tspan, trap_params, atom_params, α, N; freq=10, skip=1000)
#     samples, acc_rate = boltzmann_samples(trap_params, atom_params, N; freq=freq, skip=skip);
    
#     recapture = zeros(length(tspan));
    
#     for i ∈ 1:N
#         recapture += evolve(tspan, samples[i], atom_params, trap_params);
#     end;
    
#     return recapture ./ N, acc_rate
# end;