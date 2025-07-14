#Laser phase noise spectral density
function Sϕ(f, laser_phase_params)
    h0, hg, σg, fg = laser_phase_params;
    res = 2.0 * h0 * ones(length(f));
    
    if length(hg) > 0
        for i ∈ [1:length(hg);]
            res = res .+ 2*hg[i] .* exp.(-(f .- fg[i]).^2 ./ (2 * σg[i]^2));
        end;
    end;
        
    return res ./ (f .^ 2)
end;

#Laser phase noise amplitudes
function ϕ_amplitudes(f, laser_phase_params)
    h0, hg, σg, fg = laser_phase_params;
    Δf = f[2]-f[1];
    res = 2.0 * h0 * ones(length(f));
    
    if length(hg) > 0
        for i ∈ [1:length(hg);]
            res = res .+ 2*hg[i] .* exp.(-(f .- fg[i]).^2 ./ (2 * σg[i]^2));
        end;
    end;
    
    return 2.0 .* sqrt.(Δf * res) ./ f;
end;

#Phase noise trajectory
function ϕ(tspan, f, amplitudes)
    N = length(f);
    ϕf = rand(Uniform(0.0, 2.0*π), N); #generate random phases for components
    res = vec(sum(amplitudes .* cos.(2*π * f .* tspan' .+ ϕf), dims=1));

    return res
end;