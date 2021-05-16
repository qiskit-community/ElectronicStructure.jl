####
#### Atom
####

"""
    Atom{T}

Represents a single atom and its spatial position of type `CoordT`.

# Examples
```jldoctest
julia> Atom(:Li, (0.0, 0.0, 1.4))
Atom{Float64}(:Li, (0.0, 0.0, 1.4))
```
"""
struct Atom{CoordT}
    species::Symbol
    coords::NTuple{3, CoordT}

    function Atom(species::Symbol, coords::NTuple{3, T}) where T
        haskey(PeriodicTable.elements, species) || throw(ErrorException("No element \"$species\" exists."))
        return new{T}(species, coords)
    end
end

####
#### Geometry
####

"""
    Geometry{CoordT}

Represents the geometry of a molecule, a `Vector` of `Atom{CoordT}`s, each of
which represents an atomic species and its position. `CoordT` is the type of the
coordinates.
"""
struct Geometry{CoordT}
    atoms::Vector{Atom{CoordT}}
end

"""
    Geometry()

Create an empty `Geometry` object with `Float64` coordinates.

# Examples
```jldoctest
julia> Geometry()
Geometry{Float64}(Atom{Float64}[])
```
"""
Geometry() = Geometry(Atom{Float64}[])

"""
    Geometry(atoms::Atom...)

Initialize a `Geometry` object with `atoms`.

# Examples
```jldoctest
julia> Geometry(Atom(:H, (0., 0., 0.)), Atom(:H, (0., 0., 0.7414)))
Geometry{Float64}(Atom{Float64}[Atom{Float64}(:H, (0.0, 0.0, 0.0)), Atom{Float64}(:H, (0.0, 0.0, 0.7414))])
```
"""
Geometry(atoms::Atom...) = Geometry([atoms...])

for f in (:push!, :length, :getindex, :iterate)
    @eval Base.$f(g::Geometry, args...) = $f(g.atoms, args...)
end

####
#### Geometry
####

"""
    struct MolecularSpec{CoordT}

Contains data specifying an electronic structure problem. The results of electronic
structure calculations, for example, integrals are *not* stored here. They are stored
in the [`MolecularData`](@ref) type.

# Keywords:
- `geometry::Geometry`
-  `multiplicity::Int = 1`
-  `charge::Int = 0`
-  `basis::String = "sto-3g"`

# Additional properties
-  `spin`

# Examples
```jldoctest
julia> geometry = Geometry(Atom(:H, (0., 0., 0.)), Atom(:H, (0., 0., 0.7414)));

julia> MolecularSpec(geometry=geometry)
MolecularSpec{Float64}(Geometry{Float64}(Atom{Float64}[Atom{Float64}(:H, (0.0, 0.0, 0.0)), Atom{Float64}(:H, (0.0, 0.0, 0.7414))]), 1, 0, "sto-3g")
```
"""
Base.@kwdef struct MolecularSpec{CoordT}
    geometry::Geometry{CoordT}  #  Should we use just `Vector{Atom{CoordT}}` ?
    multiplicity::Int = 1
    charge::Int = 0
    basis::String = "sto-3g"
end

function Base.getproperty(mol::MolecularSpec, sym::Symbol)
    if sym == :spin
        return mol.multiplicity - 1  # from OpenFermion
    else
        return getfield(mol, sym)
    end
end
