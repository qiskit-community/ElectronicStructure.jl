export to_pyscf

## In order to make this a global const it is defined outside the
## try/catch block.
const pyscf = PyCall.PyNULL()

try
    copy!(pyscf, PyCall.pyimport("pyscf"))
catch
    error("Unable to load python module pyscf.")
end

"""
    to_pyscf(atom::Atom)

Convert `atom` to pyscf geometry string for a single atom (or component).
"""
to_pyscf(atom::Atom) = string(atom.species, " ", join(string.(atom.coords), " "))

"""
    to_pyscf(geom::Geometry)

Convert `geom` to pyscf geometry string, a semi-colon separated list
of specifications of species and their positions.
"""
to_pyscf(geom::Geometry) = join((to_pyscf(atom) for atom in geom), ";")

function to_pyscf(md::MolecularData)
    return pyscf.gto.Mole(atom=to_pyscf(md.geometry), basis=md.basis)
end
