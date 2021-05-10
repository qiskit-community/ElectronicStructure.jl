export to_pyscf

## In order to make this a global const it is defined outside the
## try/catch block. See the `PyCall` documentation for this idiom.
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

Convert `geom` to pyscf geometry string, a semi-colon separated list of specifications of
species and their positions.
"""
to_pyscf(geom::Geometry) = join((to_pyscf(atom) for atom in geom), ";")

"""
    to_pyscf(mol_data::MolecularData)

Convert `mol_data` to a `pyscf` python `Mol` object, and call the `Mol.build()` method.
"""
function to_pyscf(md::MolecularData)
    pymol = pyscf.gto.Mole(atom=to_pyscf(md.geometry), basis=md.basis)
    pymol.spin = md.spin
    pymol.charge = md.charge
    pymol.symmetry = false
    pymol.build()
    return pymol
end
