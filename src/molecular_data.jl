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

## Remember `deepcopy` uses a good fallback. This is only shallow
function Base.copy(md::MolecularData)
    return MolecularData(md.spec, md.nuclear_repulsion,
                  md.one_body_integrals,
                  md.two_body_integrals)
end

function Base.:(==)(md1::MolecularData, md2::MolecularData)
    return md1.spec == md2.spec &&
        md1.nuclear_repulsion == md2.nuclear_repulsion &&
        md1.one_body_integrals == md2.one_body_integrals &&
        md1.two_body_integrals == md2.two_body_integrals
end

function ZChop.zchop!(md::MolecularData)
    ZChop.zchop!(md.one_body_integrals)
    ZChop.zchop!(md.two_body_integrals)
    return md
end

function ZChop.zchop(md::MolecularData)
    return MolecularData(md.spec,
                         ZChop.zchop(md.nuclear_repulsion),
                         ZChop.zchop(md.one_body_integrals),
                         ZChop.zchop(md.two_body_integrals))
end
