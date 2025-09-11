module PBWDeformations

using Preferences: Preferences, @load_preference, @set_preferences!

import ProgressMeter

SHOW_PROGRESS_BARS() = parse(Bool, get(ENV, "SHOW_PROGRESS_BARS", "true"))

using Oscar
using Oscar.AbstractAlgebra.PrettyPrinting

using Oscar.AbstractAlgebra: ProductIterator
using Oscar.AbstractAlgebra: WeakKeyIdDict

using Oscar: IntegerUnion
using Oscar: GSetByElements

using Oscar: _is_dual
using Oscar: _is_direct_sum
using Oscar: _is_exterior_power
using Oscar: _is_symmetric_power
using Oscar: _is_tensor_power
using Oscar: _is_tensor_product
using Oscar: _is_standard_module

using Oscar.LieAlgebras: combinations
using Oscar.LieAlgebras: multicombinations
using Oscar.LieAlgebras: permutations
using Oscar.LieAlgebras: permutations_with_sign

using Base: Fix1
using Base: Fix2


import Oscar.AbstractAlgebra: change_base_ring
import Oscar.AbstractAlgebra: coefficient_ring
import Oscar.AbstractAlgebra: elem_type
import Oscar.AbstractAlgebra: gen
import Oscar.AbstractAlgebra: gens
import Oscar.AbstractAlgebra: ngens
import Oscar.AbstractAlgebra: parent_type

import Oscar: base_lie_algebra
import Oscar: comm
import Oscar: data
import Oscar: edges
import Oscar: gset_by_type
import Oscar: n_edges
import Oscar: neighbors
import Oscar: n_vertices
import Oscar: simplify
import Oscar: vertices

@static if !isdefined(Oscar, :Serialization) # introduced in https://github.com/oscar-system/Oscar.jl/pull/5024 in Oscar v1.5
    Oscar.@import_all_serialization_functions
else
    using Oscar.Serialization
    import Oscar.Serialization: load_object, save_object, type_params
end

import Base: deepcopy_internal
import Base: hash
import Base: isone
import Base: iszero
import Base: length
import Base: one
import Base: parent
import Base: show
import Base: sum
import Base: zero


export AbstractLieAlgebra, AbstractLieAlgebraElem
export ArcDiagBasedDeformBasis
export ArcDiagDeformBasis
export ArcDiagram
export ArcDiagramDirected
export ArcDiagramUndirected
export ArcDiagramVertex
export DeformationMap
export DeformBasis
export GlnGraph
export GlnGraphDeformBasis
export LieAlgebra, LieAlgebraElem
export LieAlgebraModule, LieAlgebraModuleElem
export LinearCombination
export LinearLieAlgebra, LinearLieAlgebraElem
export PseudographDeformBasis
export PseudographLabelled
export SmashProductLie, SmashProductLieElem
export SmashProductLieDeform, SmashProductLieDeformElem
export StdDeformBasis

export abstract_module
export acting_group_with_sgn
export all_arc_diagrams
export all_pbwdeformations
export all_pseudographs
export arc_diagram
export base_lie_algebra
export base_module
export deform
export edge_labels
export edges
export exterior_power
export exterior_power_obj
export general_linear_lie_algebra
export get_show_colorful_html
export inneighbor
export inneighbors
export is_crossing_free
export is_in_span
export is_in_span_with_relation
export is_linearly_independent
export is_linearly_independent_with_relations
export is_pbwdeformation
export is_span_equal
export is_span_subset
export isomorphic_module_with_simple_structure
export linear_combination
export lookup_params
export lower_vertex, is_lower_vertex
export lower_vertices
export matrix_repr_basis
export n_edges
export n_lower_vertices
export n_upper_vertices
export n_vertices
export neighbor
export neighbors
export outneighbor
export outneighbors
export pbwdeform_eqs
export set_show_colorful_html
export simple_module
export smash_product
export special_linear_lie_algebra
export special_orthogonal_lie_algebra
export standard_module
export symmetric_deformation
export symmetric_power
export symmetric_power_obj
export symmetrize
export tensor_power
export tensor_power_obj
export underlying_algebra
export upper_vertex, is_upper_vertex
export upper_vertices
export vertex_index
export vertices

const VERSION_NUMBER = Base.pkgversion(@__MODULE__)

function __init__()
    if ccall(:jl_generating_output, Cint, ()) == 0 # use !Base.generating_output() from Julia 1.11
        # modifies Oscar global variables, so must be called from the __init__ function
        patch_oscar_serialization_namespace()
        register_serialization_types()
    end

    add_verbosity_scope(:PBWDeformations)
    add_verbosity_scope(:PBWDeformationsDatabase)
end

const SO = Val{:special_orthogonal}
const GL = Val{:general_linear}

include("Types.jl")

include("OscarPatches.jl")

include("LinearIndependence.jl")
include("LinearCombination.jl")
include("ModuleSimpleStructure.jl")
include("Misc.jl")

include("DeformationBases/DeformBasis.jl")

include("SmashProductLie.jl")
include("SmashProductLieDeform.jl")
include("SmashProductPBWDeformLie.jl")
include("ArcDiagram.jl")
include("Pseudograph.jl")
include("GlnGraph.jl")

include("DeformationBases/ActingGroup.jl")
include("DeformationBases/ArcDiagBasedDeformBasis.jl")
include("DeformationBases/ArcDiagDeformBasis.jl")
include("DeformationBases/PseudographDeformBasis.jl")
include("DeformationBases/GlnGraphDeformBasis.jl")
include("DeformationBases/StdDeformBasis.jl")

include("Serialization.jl")
include("Database.jl")


###############################################################################
#
# The following function stubs' actual implementations are in the folder `ext/TestExt/`.
#
###############################################################################

function test_save_load_roundtrip end

end # module
