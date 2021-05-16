"""
    to_openfermion(atom::Atom)

Convert `atom` to OpenFermion format.

# Examples
```jldoctest
julia> to_openfermion(Atom(:Li, (0.0, 0.0, 1.4)))
(:Li, (0.0, 0.0, 1.4))
```
"""
function to_openfermion(atom::Atom)
    return (atom.species, atom.coords)
end

"""
    to_openfermion(geom::Geometry)

Convert `geom` to OpenFermion format.

# Examples
```jldoctest
julia> to_openfermion(Geometry(Atom(:H, (0., 0., 0.)), Atom(:H, (0., 0., 0.7414))))
2-element Vector{Any}:
 (:H, (0.0, 0.0, 0.0))
 (:H, (0.0, 0.0, 0.7414))
```
"""
function to_openfermion(geom::Geometry)
    return Any[to_openfermion(atom) for atom in geom]  # This need not be `Any`, but we do this for clarity
end
