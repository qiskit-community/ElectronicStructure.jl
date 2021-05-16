using ElectronicStructure

using ElectronicStructure: one_electron_integrals, two_elecron_integrals,
   MolecularData, PySCF

@pyimport qiskit_nature.drivers as drivers
@pyimport openfermion
@pyimport openfermionpyscf

## using PyCall will trigger loading pyscf-specific code

using PyCall

geoms = (
    Geometry(Atom(:H, (0., 0., 0.)), Atom(:H, (0., 0., 0.7414))),

    Geometry(Atom(:Li, (0., 0., 0.)), Atom(:H, (0., 0., 1.4))),

    Geometry(Atom(:O, (0., 0., 0.)), Atom(:H, (0.757, 0.586, 0.)),
             Atom(:H, (-0.757, 0.586, 0.)))
)

## Choose one of the geometries
geom = geoms[1]
basis = "sto-3g"

"""
    get_pyscf_molecule(geom, basis)

Use the `ElectronicStructure` interface to pyscf, which tries to do as little transformation
as possible. That is, this is a direct interface to pyscf.
"""
function get_pyscf_molecule(geom, basis)
    ## Construct specification of electronic structure problem
    mol_spec = MolecularSpec(geometry=geom)

    ## Do calculations and populate MolecularData with results
    mol_pyscf = MolecularData(PySCF, mol_spec)

    return mol_pyscf
end

mol_pyscf = get_pyscf_molecule(geom, basis);

function get_openfermion_molecule(geom, basis)
    multiplicity = 1
    mol_openfermion = openfermion.chem.MolecularData(to_openfermion(geom), basis, multiplicity)
    mol_openfermion = openfermionpyscf.run_pyscf(mol_openfermion)
    return mol_openfermion
end

mol_openfermion = get_openfermion_molecule(geom, basis);

function get_qiskit_nature_molecule(geom, basis)
    driver = drivers.PySCFDriver(atom=to_pyscf(geom),
                                 unit=drivers.UnitsType.ANGSTROM,
                                 charge=0,
                                 spin=0,
                                 basis=basis)
    mol_nature = driver.run()
    return mol_nature
end

mol_nature = get_qiskit_nature_molecule(geom, basis);
