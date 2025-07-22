#Constants and scales
#------------------------------------
const c = ustrip(u"m/s", c_0);  #Speed of light
const kB = ustrip(u"J/K", k_B)  #Boltzmann constant
const mu = ustrip(u"kg", m_u);  #Unit of atomic mass

const E0 = kB * 1e-6;       #Characteristic energy in μK
const g0 = 9.81 * 1e-6;     #Gravity free fall acceleration
const vconst = sqrt(E0/mu); #Useful constant for kinetic energy
const r0 = 1e-6;            #Characteristic distance in m
#------------------------------------



#Potential energy in gaussian beam
function Π(cord, trap_params)
    U0, w0, z0 = trap_params;
    x, y, z, vx, vy, vz = cord;
    return U0 .* (1.0 .- A(x, y, z, w0, z0) .^2);
end;


#Potential energy in gaussian beam in harmonic approximation
function Π_Harmonic(cord, trap_params)
    U0, w0, z0 = trap_params;
    x, y, z, vx, vy, vz = cord;
    
    r2 = x .^2 .+ y .^2;
    return U0 .* (2.0*r2 ./w0^2 + (z ./z0).^2);
end;



#Kinetic energy
function K(cord, trap_params, m)
    U0, w0, z0 = trap_params;
    x, y, z, vx, vy, vz = cord;
    return m/vconst^2 *(vx .^2 + vy .^2 + vz .^2) / 2.0
end;



#Total energy
function H(cord, trap_params, m; harmonic=false)
    if harmonic
        return Π_Harmonic(cord, trap_params) .+ K(cord, trap_params, m)
    else
        return Π(cord, trap_params) .+ K(cord, trap_params, m)
    end;
end;


#Target distribution for Monte-Carlo
function prob_boltzmann(cord, trap_params, atom_params; harmonic=false)    
    m, T = atom_params;
    return exp.(- H(cord, trap_params, m; harmonic) ./ T)
end;


"""
    samples_generate(trap_params, atom_params, N; freq=10, skip=1000, harmonic=true)

Generate Monte-Carlo samples of initial atom coordinates and velocities

### Input

- `trap_params` -- vector [trap depth ``U_{0}`` in ``\\mu K``, beam waist radius in ``\\mu m``, beam Rayleigh length in ``\\mu m``]
- `atom_params` -- vector [atom mass in a.u., atom temperature in ``\\mu K``]
- `N` -- number of Monte-Carlo samples, the same as number of atoms
- `freq` -- (optional, default: `10`) number of Metropolis steps skipped between samples to reduce sample dependency
- `skip` -- (optional, default: `1000`) number of Metropolis steps skipped before the Markov Chain is considered to reach stationary distribution
- `harmonic` -- (optional, default: `true`) uses harmonic approximation of gaussian beam if set to `true`, otherwise uses Metropolis sampler

### Output

Vector of Monte-Carlo samples ``[x_{i}, y_{i}, z_{i}, vx_{i}, vy_{i}, vz_{i}]`` and acceptance rate of Metropolis sampler.
If `harmonic` set to `true`, acceptance rate is set to `1.0`

"""
function samples_generate(trap_params, atom_params, N; freq=10, skip=1000, harmonic=true, eps=1e-2)
    U0, w0, z0 = trap_params;
    m, T = atom_params;

    if harmonic
        mean = zeros(6);
        cov = Diagonal(([T*w0^2/(4.0*U0),T*w0^2/(4.0*U0),T*z0^2/(2.0*U0),vconst^2*T/m,vconst^2*T/m,vconst^2*T/m]));
        d = MvNormal(mean, cov);
        
        samples = Vector{Vector{Float64}}();
        while length(samples) < N
            cord = rand(d)  
            if H(cord, trap_params, m) < U0 * (1-eps)
                push!(samples, cord)
            end
        end
        return samples, 1
    else
        mean = zeros(6);
        vstep = vconst*sqrt(T/m);
        rstep = sqrt(T/U0)/2;
        cov = Diagonal(([w0*rstep, w0*rstep, z0*sqrt(2)*rstep, vstep, vstep, vstep]) .^ 2);
        d = MvNormal(mean, cov);
        samples = [[0.0, 0.0, 0.0, vstep/sqrt(3), vstep/sqrt(3), vstep/sqrt(3)]];
        u_acc = rand(Uniform(0.0, 1.0), N*freq + skip);
        acc_rate = 0;
        for i ∈ 1:N*freq + skip - 1
            cord_last = samples[end];
            cord_new = cord_last + rand(d);
            p_acc = prob_boltzmann(cord_new, trap_params, atom_params; harmonic)/prob_boltzmann(cord_last, trap_params, atom_params; harmonic);
            
            if p_acc > u_acc[i] && H(cord_new, trap_params, m; harmonic) < U0 * (1-eps)
                push!(samples, cord_new);
                acc_rate += 1; 
            else
                push!(samples, cord_last);
            end;
        end;

        return samples[1+skip:freq:end], acc_rate/(N*freq + skip)
    end;
end;


#Generate coordinate trajectory from Monte-Carlo initial conditions
function R(t, ri, vi, ω; free=false)
    # if free
    #     return ri + vi * t;
    # else
    #     return ri * cos(ω * t) + vi/ω * sin(ω * t);
    # end;
    return free ? ri + vi * t : ri * cos(ω * t) + vi/ω * sin(ω * t);
end;    


#Generate velocity trajectory from Monte-Carlo initial conditions
function V(t, ri, vi, ω; free=false)
    # if free
    #     return vi;
    # else
    #     return vi * cos(ω * t) - ri * ω * sin(ω * t);
    # end;
    return free ? vi : vi * cos(ω * t) - ri * ω * sin(ω * t);
end;       


"""
    get_trap_params(ωr, ωz, U0, λ; m = 87.0, dif=true)

Get trap parameters given trap frequencies, depth and wavelength. 

Trap is assumed to be formed by gaussian beam with equal radial oscillation frequencies.

### Input

- `ωr` -- radial oscillation frequency in `MHz`
- `ωr` -- longitudinal oscillation frequency in `MHz`
- `U0` -- trap depth in `μK`
- `λ` -- trap wavelength in `μm`
- `m` -- atom mass in atomic units
- `dif` -- (optional, default: `true`) if set to `true`, calculates trap parameters without using `U0` and outputs trap depth
### Output

Vector of trap parameters [`w0`, `z0`, `U`]

"""
function get_trap_params(ωr, ωz, U0, λ; m = 87.0, dif=true)
    if dif
        Δω = ωr - ωz;
        w0 = λ/(1.0 - Δω/ωr) /(sqrt(2.0) * π);
        z0 = π*w0^2 / λ;
        Utemp = m * (ωr / (2.0 * vconst))^2;
    else
        w0 = vconst * sqrt(4.0 * U0 / m) / ωr ;
        z0 = vconst * sqrt(2.0 * U0 / m) / ωz ;
        Utemp = U0;
    end;
    return w0, z0, Utemp
end;
