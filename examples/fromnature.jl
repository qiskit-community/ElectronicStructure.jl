## Get a QMolecule object starting with the pyscf driver,
## all within qiskit_nature

@pyimport qiskit_nature.drivers as drivers

function get_molecule_nature()
   # atom = "H .0 .0 .0; H .0 .0 0.7414"
    #    atom = "Li .0 .0 .0; H .0 .0 1.4"
    atom = "O 0. 0. 0.; H 0.757 0.586 0.; H -0.757 0.586 0."
    driver = drivers.PySCFDriver(atom=atom,
                         unit=drivers.UnitsType.ANGSTROM,
                         charge=0,
                         spin=0,
                         basis="sto3g")
    molecule = driver.run()
    return molecule
end

nat_mol = get_molecule_nature()
