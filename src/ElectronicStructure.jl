"""
    module ElectronicStructure

Provides an interface to electronic structure calculation packages.  The only
package supported is `pyscf`.

Note that doc strings for Python objects are available at the Julia REPL. For example

    help?> ElectronicStructure.pyscf.gto.Mole
"""
module ElectronicStructure

import PeriodicTable
import Requires
import ZChop

export Atom, Geometry, MolecularSpec

function __init__()
    Requires.@require PyCall="438e738f-606a-5dbb-bf0a-cddfbfd45ab0" begin
        include("pyscf.jl")
    end
end

include("molecular_spec.jl")
include("molecular_data.jl")

end # module ElectronicStructure
