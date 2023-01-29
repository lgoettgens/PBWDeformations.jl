module PBWDeformations

using Combinatorics
using Oscar
using SparseArrays

import AbstractAlgebra:
    AllParts,
    FreeAssAlgebra,
    FreeAssAlgElem,
    Generic,
    Partition,
    Ring,
    RingElement,
    canonical_unit,
    change_base_ring,
    gens,
    ngens

import Oscar: comm, normal_form

import Base: hash, length, show, sum, ==

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
    to_arcdiag


GAP = Oscar.GAP


include("Util.jl")
include("LieAlgebraStructConsts.jl")
include("QuadraticQuoAlgebra.jl")

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
