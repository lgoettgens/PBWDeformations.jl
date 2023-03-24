module PBWDeformations

using Combinatorics
using Oscar
using SparseArrays

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

import AbstractAlgebra.Generic: _matrix, matrix_repr, rels

import Oscar: action, comm, exterior_power, normal_form, symmetric_power

import Base: deepcopy_internal, hash, isequal, iszero, length, parent, show, sum, zero, +, -, *, ==

export AbstractLieAlgebra, AbstractLieAlgebraElem
export ArcDiagDeformBasis
export ArcDiagram
export DeformationMap
export DeformBasis
export LieAlgebra, LieAlgebraElem
export LieAlgebraAbstractModule, LieAlgebraAbstractModuleElem
export LieAlgebraExteriorPowerModule, LieAlgebraExteriorPowerModuleElem
export LieAlgebraModule, LieAlgebraModuleElem
export LieAlgebraStdModule, LieAlgebraStdModuleElem
export LieAlgebraSymmetricPowerModule, LieAlgebraSymmetricPowerModuleElem
export LieAlgebraTensorPowerModule, LieAlgebraTensorPowerModuleElem
export LinearLieAlgebra, LinearLieAlgebraElem
export Pseudograph2
export PseudographDeformBasis
export QuadraticRelations
export SmashProductDeformLie
export SmashProductLie
export StdDeformBasis

export abstract_module
export all_arc_diagrams
export all_pbwdeformations
export all_pseudographs
export base_liealgebra
export bracket
export coefficient_vector
export deform
export exterior_power
export general_linear_liealgebra
export highest_weight_module
export is_crossing_free
export is_pbwdeformation
export liealgebra
export lookup_data
export matrix_repr_basis
export nedges
export pbwdeform_eqs
export smash_product
export smash_product_lie_highest_weight
export special_linear_liealgebra
export special_orthogonal_liealgebra
export standard_module
export symmetric_deformation
export symmetric_power
export tensor_power
export to_arcdiag


GAP = Oscar.GAP


include("Util.jl")

include("LieAlgebras/LieAlgebra.jl")
include("LieAlgebras/AbstractLieAlgebra.jl")
include("LieAlgebras/LinearLieAlgebra.jl")
include("LieAlgebraModules/LieAlgebraModule.jl")
include("LieAlgebraModules/LieAlgebraAbstractModule.jl")
include("LieAlgebraModules/LieAlgebraExteriorPowerModule.jl")
include("LieAlgebraModules/LieAlgebraStdModule.jl")
include("LieAlgebraModules/LieAlgebraSymmetricPowerModule.jl")
include("LieAlgebraModules/LieAlgebraTensorPowerModule.jl")
include("GapWrapper.jl")

include("FreeAssAlgQuadraticRelations.jl")
include("DeformationBases/DeformBasis.jl")

include("SmashProductLie.jl")
include("SmashProductDeformLie.jl")
include("SmashProductPBWDeformLie.jl")
include("ArcDiagram.jl")
include("Pseudograph.jl")

include("DeformationBases/ArcDiagDeformBasis.jl")
include("DeformationBases/PseudographDeformBasis.jl")
include("DeformationBases/StdDeformBasis.jl")


end
