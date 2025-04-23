#Basis states
number_of_atoms = 2
const basis4 = NLevelBasis(5);
const g1 = nlevelstate(basis4, 1);
const p1 = nlevelstate(basis4, 2);
const r1 = nlevelstate(basis4, 3);
const gt1 = nlevelstate(basis4, 4);
const zero1 = nlevelstate(basis4, 5);

const Id = g1 ⊗ dagger(g1) .+ p1 ⊗ dagger(p1) .+ r1 ⊗ dagger(r1) .+ gt1 ⊗ dagger(gt1) .+ zero1 ⊗ dagger(zero1)

#1 Operator
const σgp1 = g1 ⊗ dagger(p1) ⊗ Id;
const σpg1 = p1 ⊗ dagger(g1) ⊗ Id;
const σpr1 = p1 ⊗ dagger(r1) ⊗ Id;
const σrp1 = r1 ⊗ dagger(p1) ⊗ Id;
const np1 = p1 ⊗ dagger(p1) ⊗ Id; 
const nr1 = r1 ⊗ dagger(r1) ⊗ Id;
const zz1 = zero1 ⊗ dagger(zero1) ⊗ Id ;
#decay
const σgtp1 = gt1 ⊗ dagger(p1) ⊗ Id;
const σpgt1 = p1 ⊗ dagger(gt1) ⊗ Id;

#2 Operators
const σgp2 = Id ⊗ g1 ⊗ dagger(p1);
const σpg2 = Id ⊗ p1 ⊗ dagger(g1);
const σpr2 = Id ⊗ p1 ⊗ dagger(r1);
const σrp2 = Id ⊗ r1 ⊗ dagger(p1);
const np2 = Id ⊗ p1 ⊗ dagger(p1);
const nr2 = Id ⊗ r1 ⊗ dagger(r1);
#decay
const σgtp2 = Id ⊗ gt1 ⊗ dagger(p1);
const σpgt2 = Id ⊗ p1 ⊗ dagger(gt1);
const rr = r1 ⊗ dagger(r1) ⊗ r1 ⊗ dagger(r1);
const zz2 = Id ⊗ zero1 ⊗ dagger(zero1);
#direct 
const σgr1 = g1 ⊗ dagger(r1) ⊗ Id;
const σrg1 = r1 ⊗ dagger(g1) ⊗ Id;
const σgr2 = Id ⊗ g1 ⊗ dagger(r1);
const σrg2 = Id ⊗ r1 ⊗ dagger(g1);

const two_atom_operators = [zz1, zz2, np1, nr1, σgp1, σpg1, σpr1, σrp1, np2, nr2, σgp2, σpg2, σpr2, σrp2, rr];
const direct_operators = [zz1, zz2, nr1, nr2, σgr1, σrg1, σgr2, σrg2, rr]

function phase1(t)
    if t > 0.68 #0.2222
        return 3.902
    else
        return 0
    end;
end;

#Two-photon Rydberg hamiltonian for 1 atom
"""function Hamiltonian2(Ωr1, Ωb1, Δ1, δ1, Ωr2, Ωb2, Δ2, δ2, V_rr, d)
    return TimeDependentSum(
        [
            t -> d,
            t -> d,
            t -> -Δ1(t),
            t -> -δ1(t),
            t -> Ωr1(t) ./2.0,
            t -> conj.(Ωr1(t)) ./2.0,
            t -> Ωb1(t)/2.0,
            t -> conj.(Ωb1(t)) ./2.0,

            t -> -Δ2(t),
            t -> -δ2(t),
            t -> Ωr2(t) ./2.0,
            t -> conj.(Ωr2(t)) ./2.0,
            t -> Ωb2(t)/2.0,
            t -> conj.(Ωb2(t)) ./2.0,
            t -> V_rr
        ],
        [
        np1,
        nr1, 
        σgp1, 
        σpg1, 
        σpr1, 
        σrp1, 
        np2, 
        nr2, 
        σgp2, 
        σpg2, 
        σpr2, 
        σrp2, 
        rr    
        ];
    )
end;
#Jump operators for master equation 
function JumpOperators2(decay_params)
    Γg, Γgt = decay_params;
    return [sqrt(Γg)*σgp1, sqrt(Γgt)*σgtp1, sqrt(Γg)*σgp2, sqrt(Γgt)*σgtp2], [sqrt(Γg)*σpg1, sqrt(Γgt)*σpgt1,sqrt(Γg)*σgp2, sqrt(Γgt)*σgtp2]
end;"""

function two_atom_simulation(tspan, ψ0, 
    atom_params,    trap_params,     
    samples_first, samples_second, 
    red_laser_params_first,    blue_laser_params_first,
    red_laser_params_second,    blue_laser_params_second,
    detuning_params_first,     detuning_params_second,     decay_params, V_rr;
    atom_motion=true,     free_motion=true,
    spontaneous_decay=true,     parallel=false
    )
    
    d = 0 #6800 #MHz
    N = length(samples_first); 
    if N != length(samples_second)
        println("lengths are not equal")
    end;
    ωr, ωz = trap_frequencies(atom_params, trap_params);
    Δ01, δ01 = detuning_params_first;
    Δ02, δ02 = detuning_params_second;

    if spontaneous_decay
        decay_params_temp = decay_params;
    else
        decay_params_temp = [0.0, 0.0];
    end;
    Γg, Γgt = decay_params_temp;
    J = [sqrt(Γg)*σgp1, sqrt(Γgt)*σgtp1, sqrt(Γg)*σgp2, sqrt(Γgt)*σgtp2];
    Jdagger = [sqrt(Γg)*σpg1, sqrt(Γgt)*σpgt1,sqrt(Γg)*σgp2, sqrt(Γgt)*σgtp2];

    #J, Jdagger = [sqrt(Γg)*σgp, sqrt(Γgt)*σgtp], [sqrt(Γg)*σpg, sqrt(Γgt)*σpgt];

    ρ0 = ψ0 ⊗ dagger(ψ0);

    #Density matrix averaged over realizations of laser noise and atom dynamics.
    ρ_mean = [zero(ψ0 ⊗ dagger(ψ0)) for _ ∈ 1:length(tspan)];
    ρ_temp = [zero(ψ0 ⊗ dagger(ψ0)) for _ ∈ 1:length(tspan)];

    #Second moment for error estimation of level populations. 
    ρ2_mean = [zero(ψ0 ⊗ dagger(ψ0)) for _ ∈ 1:length(tspan)];
    
    for i ∈ 1:N
        if atom_motion
            #Atom initial conditions
            x1i, y1i, z1i, vx1i, vy1i, vz1i = samples_first[i];
            x2i, y2i, z2i, vx2i, vy2i, vz2i = samples_second[i];
        else
            x1i, y1i, z1i, vx1i, vy1i, vz1i = zeros(6);
            x2i, y2i, z2i, vx2i, vy2i, vz2i = zeros(6);
        end;
        
        #Atom trajectories
        X1 = t -> R(t, x1i, vx1i, ωr; free=free_motion);
        Y1 = t -> R(t, y1i, vy1i, ωr; free=free_motion);
        Z1 = t -> R(t, z1i, vz1i, ωz; free=free_motion);
        Vz1 = t -> V(t, z1i, vz1i, ωz; free=free_motion);
        X2 = t -> R(t, x2i, vx2i, ωr; free=free_motion);
        Y2 = t -> R(t, y2i, vy2i, ωr; free=free_motion);
        Z2 = t -> R(t, z2i, vz2i, ωz; free=free_motion);
        Vz2 = t -> V(t, z2i, vz2i, ωz; free=free_motion);
        #phase = t -> eps(t);

        #Hamiltonian params trajectories
        Ht = TimeDependentSum(
        [
            t -> d,
            t -> d,
            t -> -Δ(Vz1(t), red_laser_params_first) - Δ01,
            t -> -δ(Vz1(t), red_laser_params_first, blue_laser_params_first; parallel=parallel) - δ01,
            t -> Ω(X1(t), Y1(t), Z1(t), red_laser_params_first ) / 2.0,
            t -> conj(Ω(X1(t), Y1(t), Z1(t), red_laser_params_first ) / 2.0),
            t -> exp(1.0im * eps(t)) * Ω(X1(t), Y1(t), Z1(t), blue_laser_params_first ) / 2.0,
            t -> conj(exp(1.0im * eps(t)) * Ω(X1(t), Y1(t), Z1(t), blue_laser_params_first ) / 2.0),

            t -> -Δ(Vz2(t), red_laser_params_first) - Δ02,
            t -> -δ(Vz2(t), red_laser_params_first, blue_laser_params_first; parallel=parallel) - δ02,
            t -> Ω(X2(t), Y2(t), Z2(t), red_laser_params_second ) / 2.0,
            t -> conj(Ω(X2(t), Y2(t), Z2(t), red_laser_params_second ) / 2.0),
            t -> exp(1.0im * eps(t)) * Ω(X2(t), Y2(t), Z2(t), blue_laser_params_second ) / 2.0,
            t -> conj(exp(1.0im * eps(t)) * Ω(X2(t), Y2(t), Z2(t), blue_laser_params_second ) / 2.0),

            t -> V_rr
        ],
        two_atom_operators
        );

        #Returns hamiltonian and jump operators in a form required by timeevolution.master_dynamic
        function super_operator(t, rho)
            return Ht, J, Jdagger;
        end;
        
        _, ρ_temp = timeevolution.master_dynamic(tspan, ρ0, super_operator);

        ρ_mean = ρ_mean + ρ_temp;
        ρ2_mean = ρ2_mean + ρ_temp .^ 2;
    end;

    return ρ_mean/N, ρ2_mean/N
end;


function direct_CZ_simulation(tspan, ψ0, 
    atom_params,    trap_params,     
    samples_first, samples_second, 
    red_laser_params_first,    blue_laser_params_first,
    red_laser_params_second,    blue_laser_params_second,
    detuning_params_first,     detuning_params_second,     V_rr;
    atom_motion=true,     free_motion=true,     parallel=false
    )
    
    d = 0 #6800 #MHz
    N = length(samples_first); 
    if N != length(samples_second)
        println("lengths are not equal")
    end;
    ωr, ωz = trap_frequencies(atom_params, trap_params);
    Δ01, δ01 = detuning_params_first;
    Δ02, δ02 = detuning_params_second;
    Ω0 = blue_laser_params_first[1];

    if false
        decay_params_temp = decay_params;
    else
        decay_params_temp = [0.0, 0.0];
    end;
    Γg, Γgt = decay_params_temp;
    J = [sqrt(Γg)*σgp1, sqrt(Γgt)*σgtp1, sqrt(Γg)*σgp2, sqrt(Γgt)*σgtp2];
    Jdagger = [sqrt(Γg)*σpg1, sqrt(Γgt)*σpgt1,sqrt(Γg)*σgp2, sqrt(Γgt)*σgtp2];

    #J, Jdagger = [sqrt(Γg)*σgp, sqrt(Γgt)*σgtp], [sqrt(Γg)*σpg, sqrt(Γgt)*σpgt];

    ρ0 = ψ0 ⊗ dagger(ψ0);

    #Density matrix averaged over realizations of laser noise and atom dynamics.
    ρ_mean = [zero(ψ0 ⊗ dagger(ψ0)) for _ ∈ 1:length(tspan)];
    ρ_temp = [zero(ψ0 ⊗ dagger(ψ0)) for _ ∈ 1:length(tspan)];

    #Second moment for error estimation of level populations. 
    ρ2_mean = [zero(ψ0 ⊗ dagger(ψ0)) for _ ∈ 1:length(tspan)];
    
    for i ∈ 1:N
        if atom_motion
            #Atom initial conditions
            x1i, y1i, z1i, vx1i, vy1i, vz1i = samples_first[i];
            x2i, y2i, z2i, vx2i, vy2i, vz2i = samples_second[i];
        else
            x1i, y1i, z1i, vx1i, vy1i, vz1i = zeros(6);
            x2i, y2i, z2i, vx2i, vy2i, vz2i = zeros(6);
        end;
        
        #Atom trajectories
        X1 = t -> R(t, x1i, vx1i, ωr; free=free_motion);
        Y1 = t -> R(t, y1i, vy1i, ωr; free=free_motion);
        Z1 = t -> R(t, z1i, vz1i, ωz; free=free_motion);
        Vz1 = t -> V(t, z1i, vz1i, ωz; free=free_motion);
        X2 = t -> R(t, x2i, vx2i, ωr; free=free_motion);
        Y2 = t -> R(t, y2i, vy2i, ωr; free=free_motion);
        Z2 = t -> R(t, z2i, vz2i, ωz; free=free_motion);
        Vz2 = t -> V(t, z2i, vz2i, ωz; free=free_motion);
        #phase = t -> eps(t);

        #Hamiltonian params trajectories
        Ht = TimeDependentSum(
        [
            t -> d,
            t -> d,
            t -> Δ01,
            t -> Δ02,
            
            t -> exp(1.0im * phase1(t)) * Ω0 /2.0 , #Ω(X1(t), Y1(t), Z1(t), blue_laser_params_first ) / 2.0,
            t -> conj(exp(1.0im * phase1(t)) * Ω0 /2.0),  #Ω(X1(t), Y1(t), Z1(t), blue_laser_params_first ) / 2.0),

            t -> exp(1.0im * phase1(t)) * Ω0 /2.0 , #Ω(X2(t), Y2(t), Z2(t), blue_laser_params_second ) / 2.0,
            t -> conj(exp(1.0im * phase1(t)) * Ω0 /2.0), #Ω(X2(t), Y2(t), Z2(t), blue_laser_params_second ) / 2.0),

            t -> V_rr
        ],
        direct_operators
        #[zz1, zz2, nr1, nr2, σgr1, σrg1,  σgr2, σrg2, rr]
        );

        #Returns hamiltonian and jump operators in a form required by timeevolution.master_dynamic
        function super_operator(t, rho)
            return Ht, J, Jdagger;
        end;
        
        _, ρ_temp = timeevolution.master_dynamic(tspan, ρ0, super_operator);

        ρ_mean = ρ_mean + ρ_temp;
        ρ2_mean = ρ2_mean + ρ_temp .^ 2;
    end;

    return ρ_mean/N, ρ2_mean/N
end