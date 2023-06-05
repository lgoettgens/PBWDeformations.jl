module PBWDeformations

using Oscar

import AbstractAlgebra:
    @attributes,
    @attr,
    @enable_all_show_via_expressify,
    AllParts,
    CacheDictType,
    FreeAssAlgebra,
    FreeAssAlgElem,
    FPModule,
    FPModuleElem,
    Generic,
    MatElem,
    NCRing,
    Partition,
    ProductIterator,
    Ring,
    RingElement,
    base_ring,
    canonical_unit,
    change_base_ring,
    dim,
    elem_type,
    expressify,
    gen,
    gens,
    get_attribute,
    get_attribute!,
    get_cached!,
    has_attribute,
    ngens,
    parent_type,
    set_attribute!,
    symbols

import AbstractAlgebra.Generic: _matrix, rels

import Oscar: action, comm, exterior_power, simplify, simplify!, symmetric_power

import Base: deepcopy_internal, hash, isequal, isone, iszero, length, one, parent, show, sum, zero

import Oscar.LieAlgebras:
    AbstractLieAlgebra,
    AbstractLieAlgebraElem,
    LieAlgebra,
    LieAlgebraElem,
    LieAlgebraModule,
    LieAlgebraModuleElem,
    LinearLieAlgebra,
    LinearLieAlgebraElem,
    abstract_module,
    base_lie_algebra,
    combinations,
    exterior_power,
    general_linear_lie_algebra,
    highest_weight_module,
    is_exterior_power,
    is_standard_module,
    is_symmetric_power,
    is_tensor_power,
    lie_algebra,
    matrix_repr_basis,
    multicombinations,
    permutations,
    permutations_with_sign,
    special_linear_lie_algebra,
    special_orthogonal_lie_algebra,
    standard_module,
    symmetric_power,
    tensor_power

export AbstractLieAlgebra, AbstractLieAlgebraElem
export ArcDiagDeformBasis
export ArcDiagram
export DeformationMap
export DeformBasis
export LieAlgebra, LieAlgebraElem
export LieAlgebraModule, LieAlgebraModuleElem
export LinearLieAlgebra, LinearLieAlgebraElem
export Pseudograph2
export PseudographDeformBasis
export QuadraticRelations
export SmashProductLie, SmashProductLieElem
export SmashProductLieDeform, SmashProductLieDeformElem
export StdDeformBasis

export abstract_module
export all_arc_diagrams
export all_pbwdeformations
export all_pseudographs
export base_lie_algebra
export deform
export exterior_power
export flatten
export general_linear_lie_algebra
export highest_weight_module
export is_crossing_free
export is_pbwdeformation
export lie_algebra
export lie_module
export lookup_data
export matrix_repr_basis
export nedges
export pbwdeform_eqs
export smash_product
export special_linear_lie_algebra
export special_orthogonal_lie_algebra
export standard_module
export symmetric_deformation
export symmetric_power
export tensor_power
export to_arcdiag
export underlying_algebra

function __init__()
    Hecke.add_verbose_scope(:PBWDeformations)
end

include("Util.jl")

include("FreeAssAlgQuadraticRelations.jl")
include("DeformationBases/DeformBasis.jl")

include("SmashProductLie.jl")
include("SmashProductLieDeform.jl")
include("SmashProductPBWDeformLie.jl")
include("ArcDiagram.jl")
include("Pseudograph.jl")
include("OldStuff.jl")

include("DeformationBases/ArcDiagDeformBasis.jl")
include("DeformationBases/PseudographDeformBasis.jl")
include("DeformationBases/StdDeformBasis.jl")


end
