####
#### Atom
####

"""
    Atom{T}

Represents a single atom and it's spatial position of type `CoordT`.
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

Represents the geometry of a molecule, a `Vector` of `Atom{CoordT}`s,
each of which represents an atomic species and its position.
"""
struct Geometry{CoordT}
    atoms::Vector{Atom{CoordT}}
end

"""
    Geometry()

Create an empty `Geometry` object with `Float64` coordinates.
"""
Geometry() = Geometry{Atom{Float64}}(Vector{Atom{Float64}}[])

for f in (:push!, :length, :getindex, :iterate)
    @eval Base.$f(g::Geometry, args...) = $f(g.atoms, args...)
end

####
#### Geometry
####

"""
    struct MolecularData{CoordT}

Contains data specifying an electronic structure problem.
"""
Base.@kwdef struct MolecularData{CoordT}
    geometry::Geometry{CoordT}
    multiplicity::Int = 1
    charge::Int = 0
    basis::String = "sto-3g"
end
