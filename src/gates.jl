function get_gate(U)
    op = dense(identityoperator(ColdAtoms.basis))
    op.data[1:2, 1:2] = U
    return op
end


function project_on_qubit(ρ, n = 1)
    basis_sub = SubspaceBasis(basis, [ket_0, ket_1])
    P = reduce(⊗, [projector(basis_sub, basis) for _ = 1:n])
    return P * ρ * dagger(P)
end;


Hadamard = get_gate(Matrix{ComplexF64}([1 1; 1 -1]/sqrt(2)))
X = get_gate(Matrix{ComplexF64}([
    0 1;
    1 0
]))
Y = get_gate(Matrix{ComplexF64}([
    0 -1.0im;
    1.0im 0
]));
Z = get_gate(Matrix{ComplexF64}([
    1 0;
    0 -1
]));
RX = θ -> get_gate(Matrix{ComplexF64}([cos(θ/2) -1.0im*sin(θ/2); -1.0im*sin(θ/2) cos(θ/2)]));
RY = θ -> get_gate(Matrix{ComplexF64}([cos(θ/2) -sin(θ/2); sin(θ/2) cos(θ/2)]));
RZ = θ -> get_gate(Matrix{ComplexF64}([
    exp(-1.0im*θ/2) 0;
    0 exp(1.0im*θ/2)
]));

# CNOT = Matrix{ComplexF64}([
#     1 0 0 0;
#     0 1 0 0;
#     0 0 0 1;
#     0 0 1 0]);
# CZ = Matrix{ComplexF64}([
#     1 0 0 0;
#     0 1 0 0;
#     0 0 1 0;
#     0 0 0 -1;
# ])
