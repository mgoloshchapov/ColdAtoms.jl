module ColdAtoms

using Distributions, Random
using PhysicalConstants.CODATA2018: c_0, k_B, m_u
using Unitful
using LinearAlgebra
using QuantumOptics
using SplitApplyCombine
using Interpolations
using Polynomials
using SpecialPolynomials
using HypergeometricFunctions

using OrderedCollections
using Plots
using ProgressBars
#using OrdinaryDiffEq

export 
    w0_to_z0, trap_frequencies, E, I, get_rydberg_probs,
    release_recapture,
    samples_generate, R, V, get_trap_params, H, samples_visualise,
    Sϕ, ϕ_amplitudes, ϕ,
    simulation, Ω_twophoton, T_twophoton, δ_twophoton, Ωr_required, 
    ket_0, ket_1, ket_r, ket_p, ket_l,
    simple_flattopHG_field,simple_flattopLG_field,
    HG_coeff, simulation_blue_intens,
    g1, p1, r1, gt1, zero1,
    Id, two_atom_simulation,
    direct_CZ_simulation
        
include("utilities.jl")
include("basic_experiments.jl")
include("lasernoise_sampler.jl")
include("atom_sampler.jl")
include("rydberg_model.jl")
include("arbitrary_beams.jl")
include("rydberg_model_arb_bms.jl")
include("two_atom_model.jl")

end