using Artifacts
import JLD2
using ElectronicStructure: MolecularData

## TODO: encapsulate this somehow
"""
    get_molecular_data()

Return a `Dict` containing a few examples of PySCF electronic structure
calculations. They are retrieved in `JLD2` format using the `Artifacts`
framework. The `JLD2` file is read and the result is returned.

The data returned was generated using the function `generate_pyscf_molecules`
in *compare_formats.jl*.
"""
function get_molecular_data()
    rootpath = artifact"molecular_data"
    filename = joinpath(rootpath, "molecular_data_examples.jld2")
    data = JLD2.load(filename)
    return data
end

@testset "molecular data sets" begin
    data_top = get_molecular_data()
    @test collect(keys(data_top)) == ["header", "data"]
    data = data_top["data"]
    for mol_name in keys(data)
        mol1 = data[mol_name]
        @test isa(mol1, MolecularData)
        @test find_index_order(mol1.two_body_integrals) == :chemist
        iop = InteractionOperator(mol1)
        @test find_index_order(iop.two_body_tensor) == :physicist
        @test find_index_order(iop) == :physicist
        iop2 = InteractionOperator(mol1; index_order=:chemist)
        @test find_index_order(iop2.two_body_tensor) == :chemist
        iop2 = InteractionOperator(mol1; index_order=:intermediate)
        @test find_index_order(iop2.two_body_tensor) == :intermediate
        for spin_order in (:block, :alternating, :none)
            @test find_index_order(InteractionOperator(mol1; spin_order=spin_order).two_body_tensor) == :physicist
        end
    end
end
