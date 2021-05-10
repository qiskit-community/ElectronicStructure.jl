####
#### Transformations to our representations
####

## Caution, this differs from OpenFermion by a factor of 1/2 in the two-body integral.
## In OpenFermion, they include the factor when constructing InteractionOperator
function spin_orbital_from_spatial(one_body_integrals, two_body_integrals)
    n_side = 2 * size(one_body_integrals)[1]
    r1 = 1:2:n_side # odd indices
    r2 = 2:2:n_side # even indices

    one_body_coefficients = similar(one_body_integrals, (n_side, n_side))
    one_body_coefficients[r1, r1] .= one_body_integrals
    one_body_coefficients[r2, r2] .= one_body_integrals

    two_body_coefficients = similar(two_body_integrals, (n_side, n_side, n_side, n_side))
    for inds in ((r1, r2, r2, r1), (r2, r1, r1, r2), (r1, r1, r1, r1), (r2, r2, r2, r2))
        two_body_coefficients[inds...] .= two_body_integrals
    end
    two_body_integrals .= one(eltype(two_body_integrals)) / 2 .* two_body_integrals

    return ZChop.zchop!.((one_body_coefficients, two_body_coefficients))
end

struct InteractionOperator
    nuclear_repulsion::Float64
    one_body_tensor::Array{Float64, 2}
    two_body_tensor::Array{Float64, 4}
end

## NOTE !!!! Conversion of InteractionOperator to FermionOperator is in OF conversions.py
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function InteractionOperator(mol_data::MolecularData)
    tens1, tens2 = spin_orbital_from_spatial(mol_data.one_body_integrals, mol_data.two_body_integrals)
    return InteractionOperator(mol_data.nuclear_repulsion, tens1, tens2)
end
