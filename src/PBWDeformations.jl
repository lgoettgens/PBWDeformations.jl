module PBWDeformations

using Combinatorics
using Oscar
using SparseArrays

import AbstractAlgebra:
    AllParts,
    FreeAssAlgebra,
    FreeAssAlgElem,
    Generic,
    MatElem,
    Partition,
    ProductIterator,
    Ring,
    RingElement,
    base_ring,
    canonical_unit,
    change_base_ring,
    elem_type,
    gen,
    gens,
    ngens,
    parent_type,
    symbols

import AbstractAlgebra.PrettyPrinting: @enable_all_show_via_expressify, expressify

import Oscar: comm, exterior_power, normal_form, symmetric_power

import Base: deepcopy_internal, hash, isequal, iszero, length, parent, show, sum, zero, +, -, *, ==

export ArcDiagram,
    ArcDiagDeformBasis,
    DeformationMap,
    DeformBasis,
    Pseudograph2,
    PseudographDeformBasis,
    QuadraticRelations,
    SmashProductLie,
    SmashProductLieInfo,
    SmashProductDeformLie,
    StdDeformBasis

export all_arc_diagrams,
    all_pseudographs,
    is_crossing_free,
    is_pbwdeform,
    exterior_power,
    lookup_data,
    nedges,
    pbwdeforms_all,
    pbwdeform_eqs,
    smash_product_lie,
    smash_product_lie_highest_weight,
    smash_product_lie_so_fundamental_module,
    smash_product_lie_so_symmpowers_standard_module,
    smash_product_lie_so_extpowers_standard_module,
    smash_product_lie_sp_symmpowers_standard_module,
    smash_product_lie_sp_extpowers_standard_module,
    smash_product_deform_lie,
    smash_product_symmdeform_lie,
    symmetric_power,
    so_standard_module,
    tensor_power,
    to_arcdiag


GAP = Oscar.GAP


include("Util.jl")
include("FreeAssAlgQuadraticRelations.jl")
include("LieAlgebraStructConsts.jl")
include("SOnModules/SOnModule.jl")
include("SOnModules/SOnStdModule.jl")
include("SOnModules/SOnTensorPowerModule.jl")
include("SOnModules/SOnExteriorPowerModule.jl")
include("SOnModules/SOnSymmetricPowerModule.jl")

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
