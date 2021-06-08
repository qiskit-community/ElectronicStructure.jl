using ElectronicStructure
using Artifacts
import JLD2

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
