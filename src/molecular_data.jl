struct MolecularData{T}
    spec::MolecularSpec
    nuclear_repulsion::Float64
    one_body_integrals::Array{T, 2}
    two_body_integrals::Array{T, 4}
end
