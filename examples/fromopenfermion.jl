@pyimport openfermion
@pyimport openfermionpyscf

function get_molecule_openfermion()
    basis = "sto-3g"
    multiplicity = 1
    run_scf = 1
    run_mp2 = 1
    run_cisd = 0
    run_ccsd = 0
    run_fci = 1
    bond_length = 0.7414
    geometry = [("H", (0., 0., 0.)), ("H", (0., 0., bond_length))]
    molecule = openfermion.chem.MolecularData(
        geometry, basis, multiplicity)

    molecule = openfermionpyscf.run_pyscf(molecule)
                         # run_scf=run_scf,
                         # run_mp2=run_mp2,
                         # run_cisd=run_cisd,
                         # run_ccsd=run_ccsd,
                         # run_fci=run_fci)
    return molecule
end

mol_of = get_molecule_openfermion();
