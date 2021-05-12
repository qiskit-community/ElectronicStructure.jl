using PyCall: PyObject

export PyMol, PySCF
export to_pyscf, hartree_fock!

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
    to_pyscf(mol_data::MolecularSpec)

Convert `mol_data` to a `pyscf` python `Mol` object, and call the `Mol.build()` method.
"""
function to_pyscf(md::MolecularSpec)
    pymol = pyscf.gto.Mole(atom=to_pyscf(md.geometry), basis=md.basis)
    pymol.spin = md.spin
    pymol.charge = md.charge
    pymol.symmetry = false
    pymol.build()
    return pymol
end

"""
    struct PyMol

Wraps the Python class `pyscf.gto.Mole`.
"""
struct PyMol
    pymol::PyObject

    function PyMol(pymol::PyObject)
        class = PyCall.pytypeof(pymol).__name__
        if class != "Mole"
            throw(ErrorException("`PyMol` expecting python class `Mole`."))
        end
        return new(pymol)
    end
end

"""
    PyMol(md::MolecularSpec)

Create an instance of Python class `pyscf.gto.Mole` and wrap in `PyMol`.  This essentially
converts the `ElectronicStructure`-native `MolecularSpec` type to the corresponding `pyscf`
type.
"""
PyMol(md::MolecularSpec) = PyMol(to_pyscf(md))

"""
    struct PySCF

Wraps an `scf` calculation object from `pyscf`. The `scf` is the central object
for nearly all calculations available from `pyscf`. You must call `hartree_fock!`
on this object before running other calculations.
"""
struct PySCF
    scf::PyObject
end

function PySCF(pymol::PyMol)
    if pymol.pymol.spin != 0  # from OpenFermion
        pyscf_scf = pyscf.scf.ROHF(pymol.pymol)
    else
        pyscf_scf = pyscf.scf.RHF(pymol.pymol)
    end
    return PySCF(pyscf_scf)
end

"""
    hartree_fock!(scf::PySCF)

Perform Hartee-Fock calculation using the `scf` calculational object.  This method must be
called on `scf` before any other calculations, such as integrals may be performed.
"""
function hartree_fock!(scf::PySCF)
    verbose = scf.scf.verbose
    scf.scf.verbose = 0
    scf.scf.run()
    scf.scf.verbose = verbose
    return scf
end

## These integrals, etc. are calculated in two places in OpenFermion,
## perhaps and older and newer way of doing it.
## _pyscf_molecular_data.y seems simpler, more organized, and clear
## _run_pyscf.py is more complicated.
## The latter makes a distinction between 1d compressed and not integrals

## Both OpenFermion and qiskit.nature do mo_coeff' * hcore * mo_coeff.
## OpenFermion reshapes this into one_electron, but it really looks like a no-op to me.
## Wait,... In _pyscf_molecular_data.py, they seem to agree with me.
function one_electron_integrals(scf::PySCF)
    mo_coeff = scf.scf.mo_coeff
    hcore = scf.scf.get_hcore()
    return mo_coeff' * hcore * mo_coeff
end

function two_elecron_integrals(pymol::PyMol, scf::PySCF)
    two_elecron_compressed = pyscf.ao2mo.kernel(pymol.pymol, scf.scf.mo_coeff)
    n_orbitals = size(scf.scf.mo_coeff)[2]
    symmetry_code = 1 # No permutation symmetry
    return pyscf.ao2mo.restore(symmetry_code, two_elecron_compressed, n_orbitals)
end

## See OpenFermion molecular_data.py
## get_molecular_hamiltonian for how to convert these integrals to InteractionOperator

function MolecularData(::Type{PySCF} ,mol_spec::MolecularSpec)
    pymol = PyMol(mol_spec)  # Convert MolecularSpec to pyscf equiv.
    scf = PySCF(pymol)       # Create scf calc object.
    hartree_fock!(scf)       # Run basic, neccesary calcs.
    one_e_ints = one_electron_integrals(scf)  # Compute one and two body integrals
    two_e_ints = two_elecron_integrals(pymol, scf)
    nuclear_repulsion = pymol.pymol.energy_nuc() # Compute constant energy
    return MolecularData(mol_spec, nuclear_repulsion, one_e_ints, two_e_ints)
end
