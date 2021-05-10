"""
    module ElectronicStructure

Provides an interface to electronic structure calculation packages.  The only
package supported is `pyscf`.
"""
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

end # module ElectronicStructure
