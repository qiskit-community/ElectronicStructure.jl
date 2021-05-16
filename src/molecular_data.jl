"""
    struct MolecularData{T}

Structure describing molecular Hamiltonian. This stores the specification
of the problem and the resulting integrals and nuclear repulsion energy.
"""
struct MolecularData{T}
    spec::MolecularSpec
    nuclear_repulsion::Float64
    one_body_integrals::Array{T, 2}
    two_body_integrals::Array{T, 4}
end
