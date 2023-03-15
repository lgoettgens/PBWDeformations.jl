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

export AbstractLieAlgebra,
    AbstractLieAlgebraElem,
    ArcDiagram,
    ArcDiagDeformBasis,
    DeformationMap,
    DeformBasis,
    LieAlgebra,
    LieAlgebraElem,
    LieAlgebraModule,
    LieAlgebraModuleElem,
    LieAlgebraAbstractModule,
    LieAlgebraAbstractModuleElem,
    LieAlgebraExteriorPowerModule,
    LieAlgebraExteriorPowerModuleElem,
    LieAlgebraStdModule,
    LieAlgebraStdModuleElem,
    LieAlgebraSymmetricPowerModule,
    LieAlgebraSymmetricPowerModuleElem,
    LieAlgebraTensorPowerModule,
    LieAlgebraTensorPowerModuleElem,
    LinearLieAlgebra,
    LinearLieAlgebraElem,
    Pseudograph2,
    PseudographDeformBasis,
    QuadraticRelations,
    SmashProductLie,
    SmashProductDeformLie,
    StdDeformBasis

export abstract_module,
    all_arc_diagrams,
    all_pseudographs,
    base_liealgebra,
    bracket,
    coefficient_vector,
    is_crossing_free,
    is_pbwdeform,
    exterior_power,
    general_linear_liealgebra,
    highest_weight_module,
    liealgebra,
    lookup_data,
    matrix_repr_basis,
    nedges,
    pbwdeforms_all,
    pbwdeform_eqs,
    smash_product_lie,
    smash_product_lie_highest_weight,
    smash_product_lie_so_symmpowers_standard_module,
    smash_product_lie_so_extpowers_standard_module,
    smash_product_deform_lie,
    smash_product_symmdeform_lie,
    special_linear_liealgebra,
    special_orthogonal_liealgebra,
    symmetric_power,
    standard_module,
    tensor_power,
    to_arcdiag


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
