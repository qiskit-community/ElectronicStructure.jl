####
#### Atom
####

struct Atom{T}
    species::Symbol
    coords::NTuple{3, T}

    function Atom(species::Symbol, coords::NTuple{3, T}) where T
        haskey(PeriodicTable.elements, species) || throw(ErrorException("No element \"$species\" exists."))
        return new{T}(species, coords)
    end
end

####
#### Geometry
####

struct Geometry{T}
    atoms::Vector{T}
end

Geometry() = Geometry{Atom{Float64}}(Vector{Atom{Float64}}[])

for f in (:push!, :length)
    @eval Base.$f(g::Geometry, args...) = $f(g.atoms, args...)
end

####
#### Geometry
####

Base.@kwdef struct MolecularData{T}
    geometry::Geometry{T}
    multiplicity::Int = 1
    charge::Int = 0
    basis::String = "sto-3g"
end
