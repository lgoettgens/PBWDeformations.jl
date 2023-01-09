module PBWDeformations

using Combinatorics
using Oscar
using SparseArrays

import AbstractAlgebra:
    AllParts,
    Generic,
    NCRing,
    NCRingElem,
    Partition,
    Ring,
    RingElement,
    base_ring,
    canonical_unit,
    change_base_ring,
    check_parent,
    coeff,
    divexact,
    elem_type,
    gen,
    gens,
    is_gen,
    is_monomial,
    monomial,
    ngens,
    parent_type,
    quo,
    symbols,
    vars

import Oscar: comm, normal_form

import Base:
    Array,
    deepcopy,
    deepcopy_internal,
    hash,
    isequal,
    isone,
    iszero,
    length,
    one,
    parent,
    show,
    sum,
    xor,
    zero,
    +,
    -,
    *,
    ^,
    ==

export ArcDiagram,
    ArcDiagDeformBasis,
    DeformationMap,
    DeformBasis,
    FreeAlgebra,
    FreeAlgebraElem,
    Pseudograph2,
    PseudographDeformBasis,
    QuadraticQuoAlgebra,
    QuadraticQuoAlgebraElem,
    SmashProductLie,
    SmashProductLieInfo,
    SmashProductDeformLie,
    StdDeformBasis

export all_arc_diagrams,
    all_pseudographs,
    corresponding_arc_diagram,
    corresponding_arc_diagrams,
    free_algebra,
    is_crossing_free,
    is_pbwdeform,
    lookup_data,
    nedges,
    quadratic_quo_algebra,
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
include("Algebra.jl")
include("FreeAlgebra.jl")
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
