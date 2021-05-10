module ElectronicStructure

import PeriodicTable
import Requires

export Atom, Geometry, MolecularData

function __init__()
    Requires.@require PyCall="438e738f-606a-5dbb-bf0a-cddfbfd45ab0" begin
        include("pyscf.jl")
    end
end


include("molecular_data.jl")
#include("pyscf.jl")

end # module ElectronicStructure
