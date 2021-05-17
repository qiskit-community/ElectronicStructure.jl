using ElectronicStructure

using PyCall

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

function get_openfermion_molecule(geom, basis)
    multiplicity = 1
    mol_openfermion = openfermion.chem.MolecularData(to_openfermion(geom), basis, multiplicity)
    mol_openfermion = openfermionpyscf.run_pyscf(mol_openfermion)
    return mol_openfermion
end

function get_qiskit_nature_molecule(geom, basis)
    driver = drivers.PySCFDriver(atom=to_pyscf(geom),
                                 unit=drivers.UnitsType.ANGSTROM,
                                 charge=0,
                                 spin=0,
                                 basis=basis)
    mol_nature = driver.run()
    return mol_nature
end

mol_pyscf = get_pyscf_molecule(geom, basis);
mol_openfermion = get_openfermion_molecule(geom, basis);
mol_nature = get_qiskit_nature_molecule(geom, basis);

function boolstr(x::Bool)
    if x
        return "true"
    else
        return "*** false"
    end
end

tstpr(tst, msg) = println(boolstr(tst), ": ", msg)

function compare_calculations()
    tst = mol_nature.mo_onee_ints ≈ mol_openfermion.one_body_integrals
    tstpr(tst, "(raw) nature and openfermion one-body integrals are equal")

    tst = mol_openfermion.one_body_integrals ≈ mol_pyscf.one_body_integrals
    tstpr(tst, "openfermion and pyscf one-body integrals are equal")

    tst = mol_nature.mo_eri_ints ≈ mol_pyscf.two_body_integrals
    tstpr(tst, "(raw) nature and pyscf two-body integrals are equal")

    tst = mol_openfermion.two_body_integrals ≈ chem_to_phys(mol_pyscf.two_body_integrals)
    tstpr(tst, "openfermion and phys_to_chem(pyscf) two-body integrals are equal")

    tst = phys_to_chem(mol_openfermion.two_body_integrals) ≈ mol_pyscf.two_body_integrals
    tstpr(tst, "chem_to_phys(openfermion) and pyscf two-body integrals are equal")

    tst = check_two_body_symmetries(mol_pyscf.two_body_integrals; chemist=true)
    tstpr(tst, "pyscf two-body integrals have chemists symmetry")

    tst = check_two_body_symmetries(chem_to_phys(mol_pyscf.two_body_integrals); chemist=false)
    tstpr(tst, "pyscf two-body integrals have physicists' symmetry under `chem_to_phys`")

    tst = check_two_body_symmetries(mol_nature.mo_eri_ints; chemist=true)
    tstpr(tst, "nature (raw) two-body integrals have chemists symmetry")

    tst = check_two_body_symmetries(chem_to_phys(mol_nature.mo_eri_ints); chemist=false)
    tstpr(tst, "nature (raw) two-body integrals have physicists' symmetry under `chem_to_phys`")

    tst = check_two_body_symmetries(mol_openfermion.two_body_integrals; chemist=false)
    tstpr(tst, "openfermion two-body integrals have physicists' symmetry")

    tst = check_two_body_symmetries(phys_to_chem(mol_openfermion.two_body_integrals); chemist=true)
    tstpr(tst, "openfermion two-body integrals have chemists' symmetry under `phys_to_chem`")

    ## get spin orbitals, double the number of orbitals
    molham_openfermion = mol_openfermion.get_molecular_hamiltonian()
    tst = check_two_body_symmetries(molham_openfermion.two_body_tensor; chemist=false)
    tstpr(tst, "openfermion mol. ham. two-body tensor has physicists' symmetry.")

    ## This is buggy. Gives wrong tensorsa
    # iop = InteractionOperator(mol_pyscf; block_spin=true, transform=nothing)

    tst = check_two_body_symmetries(mol_nature.two_body_integrals; chemist=false)
    tstpr(tst, "qiskit_nature spin orb two_body_tensor has PHYSICISTs symmetry.")
end
