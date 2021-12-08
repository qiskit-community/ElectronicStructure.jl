"""
    module ElectronicStructure

Provides an interface to electronic structure calculation packages.

PySCF is supported by the package ElectronicStructurePySCF
"""
module ElectronicStructure


import PeriodicTable
import Requires
import ZChop
using Tullio: @tullio

import PkgVersion
import Dates

export Atom, Geometry, MolecularSpec, InteractionOperator, InteractionOperator!
export non_zero_elements, non_zero_elements_python, phys_to_chem, chem_to_phys,
    check_two_body_symmetries, find_index_order,
    to_index_order

export to_openfermion, to_pyscf

# function __init__()
#     Requires.@require PyCall="438e738f-606a-5dbb-bf0a-cddfbfd45ab0" begin
#         include("pyscf.jl")
#     end
# end

include("utils.jl")
include("molecular_spec.jl")
include("molecular_data.jl")
include("interaction_operator.jl")
include("openfermion.jl")

end # module ElectronicStructure
