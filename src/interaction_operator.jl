# Disable this because of broken dependencies
# using FastBroadcast: @..

####
#### Transformations to our representations
####

## Caution, we choose that this differs from OpenFermion by a factor of 1/2 in the two-body integral.
## In OpenFermion, they include the factor when constructing InteractionOperator.
## See spinorb_from_spatial in molecular_data.py.
"""
    spin_orbital_from_spatial(one_body_integrals, two_body_integrals; block_spin=false, index_order=:physicist)

Create coefficients for spin-orbitals from coefficients for spatial orbitals only.

If `block_spin` is `true` the spins are in two sectors. Otherwise, the spin variable varies
more rapidly than the spatial variable, i.e. the spins are interleaved; orbitals
corresponding to the same spatial orbital, but with differing spins are adjacent.  Note that
Qiskit uses the former ordering and OpenFermion uses the latter.

`index_order` must be one of `:physicist`, `:chemist`, or `:intermediate`, and must correspond to the
actual symmetry of `two_body_integrals` in order that the spin-orbital tensor be constructed correctly.
This symmetry may be determined via `find_index_order`.
"""
function spin_orbital_from_spatial(one_body_integrals, two_body_integrals; block_spin=false, index_order=:physicist)
    spatial_dim = size(one_body_integrals)[1]
    spin_orb_dim = 2 * spatial_dim
    if block_spin
        r1 = 1:spatial_dim  # first block of orbitals
        r2 = (spatial_dim+1):spin_orb_dim  # second block of orbitals
    else  # interleaved spin order
        r1 = 1:2:spin_orb_dim # odd indices
        r2 = 2:2:spin_orb_dim # even indices
    end
    return _spin_orbital_from_spatial(one_body_integrals, two_body_integrals, r1, r2, index_order)
end

## Create spin-orbital coefficients from spatial-orbital coeffcients. r1 and r2 are indices
## encoding whether spin variable varies faster or slower than spatial-orbital variable.
function _spin_orbital_from_spatial(one_body_integrals, two_body_integrals, r1, r2, index_order=:physicist)
    new_dim = 2 * size(one_body_integrals)[1]

    one_body_coefficients = zeros(eltype(one_body_integrals), (new_dim, new_dim))
    @. one_body_coefficients[r1, r1] = one_body_integrals
    @. one_body_coefficients[r2, r2] = one_body_integrals

    two_body_coefficients = zeros(eltype(two_body_integrals), (fill(new_dim, 4)...,))
    if index_order == :physicist
        index_set = ((r1, r2, r2, r1), (r2, r1, r1, r2), (r1, r1, r1, r1), (r2, r2, r2, r2))
    elseif index_order == :intermediate
        index_set = ((r1, r2, r1, r2), (r2, r1, r2, r1), (r1, r1, r1, r1), (r2, r2, r2, r2))
    elseif index_order == :chemist
        index_set = ((r1, r1, r2, r2), (r2, r2, r1, r1), (r1, r1, r1, r1), (r2, r2, r2, r2))
    else
        throw(ArgumentError("Unknown index order: $index_order"))
    end
    @inbounds for inds in index_set
        @. two_body_coefficients[inds...] = two_body_integrals
    end
    @. two_body_coefficients = two_body_coefficients / 2

    # Should we do zchop ?
    # return ZChop.zchop!.((one_body_coefficients, two_body_coefficients))
    return (one_body_coefficients, two_body_coefficients)
end

"""
    struct InteractionOperator

Structure to represent a molecular Hamiltonian.
"""
struct InteractionOperator
    nuclear_repulsion::Float64
    one_body_tensor::Array{Float64, 2}
    two_body_tensor::Array{Float64, 4}
end

Base.:(==)(iop1::InteractionOperator, iop2::InteractionOperator) =
    iop1.nuclear_repulsion == iop2.nuclear_repulsion &&
    iop1.one_body_tensor == iop2.one_body_tensor &&
    iop1.two_body_tensor == iop2.two_body_tensor

Base.isapprox(iop1::InteractionOperator, iop2::InteractionOperator) =
    isapprox(iop1.nuclear_repulsion, iop2.nuclear_repulsion) &&
    isapprox(iop1.one_body_tensor, iop2.one_body_tensor) &&
    isapprox(iop1.two_body_tensor, iop2.two_body_tensor)

## NOTE !!!! Conversion of InteractionOperator to FermionOperator is in OF conversions.py
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
"""
    InteractionOperator(mol_data::MolecularData; spin_order=:block, index_order=:default)

Create an InteractionOperator from `mol_data`.

`mol_data` contains integrals for spatial orbitals. `InteractionOrbitals` creates coefficients for
spin orbitals, doubling the number of orbitals.

If `spin_order` is `:block` the two-body tensor is constructed with spins in two
sectors. If `spin_order` is `:alternating`, the spin variable varies more rapidly than the spatial variable,
i.e. adjacent orbitals have differing spins. Note that Qiskit uses the former ordering and
OpenFermion uses the latter. If `spin_order` is `:none`, then the operator includes only the input
spatial oribtals; no spin orbitals are computed.

If `index_order` is `:default` then it is assumed that the two-body tensor in `mol_data` is
in chemists' index order and that the desired order for the spin-orbital tensor is
physicists' order.  Otherwise, `index_order` must be one of `:chemist`, `:physicist`, or
`:intermediate` and an attempt will be made to construct the spin-orbital tensor in the
corresponding index order.

The default values of `block_spin` and `index_order` agree with the `InteractionOperator` from OpenFermion.
"""
function InteractionOperator(mol_data::MolecularData; spin_order=:block, index_order=:default)
    interaction_operator(mol_data.nuclear_repulsion, mol_data.one_body_integrals, mol_data.two_body_integrals;
                         spin_order=spin_order, index_order=index_order)
end

function interaction_operator(nuclear_repulsion, one_body_integrals, two_body_integrals;
                              spin_order=:block, index_order=:default)
    tb = two_body_integrals
    if index_order == :default
        tens_tmp = chem_to_phys(tb)
        index_order = :physicist
    else
        tens_tmp = to_index_order(tb, index_order)
    end
    if spin_order == :none
        return InteractionOperator(nuclear_repulsion, one_body_integrals, tens_tmp)
    end
    if spin_order == :block
        block_spin = true
    else
        block_spin = false
    end
    tens1, tens2 = spin_orbital_from_spatial(one_body_integrals, tens_tmp, block_spin=block_spin, index_order=index_order)
    return InteractionOperator(nuclear_repulsion, tens1, tens2)
end

"""
    zchop!(iop::InteractionOperator, zeps=ZChop.ZEPS)

Remove terms in `iop` that are close to zero.
"""
function ZChop.zchop!(iop::InteractionOperator, zeps::Real=ZChop.ZEPS)
    ZChop.zchop!(iop.one_body_tensor, zeps)
    ZChop.zchop!(iop.two_body_tensor, zeps)
    return iop
end

function ZChop.nchop(iop::InteractionOperator)
    return InteractionOperator(ZChop.nchop(iop.nuclear_repulsion),
                               ZChop.nchop(iop.one_body_tensor),
                               ZChop.nchop(iop.two_body_tensor))
end

"""
    zchop(iop::InteractionOperator, zeps=ZChop.ZEPS)

See `zchop!(iop::InteractionOperator)`
"""
function ZChop.zchop(iop::InteractionOperator, zeps::Real=ZChop.ZEPS)
    return InteractionOperator(ZChop.zchop(iop.nuclear_repulsion, zeps),
                               ZChop.zchop(iop.one_body_tensor, zeps),
                               ZChop.zchop(iop.two_body_tensor, zeps))
end

"""
    find_index_order(iop::InteractionOperator)

Return the index-ordering convention of the two-body tensor in `iop`.
"""
find_index_order(iop::InteractionOperator) = find_index_order(iop.two_body_tensor)
