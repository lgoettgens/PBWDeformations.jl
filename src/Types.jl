################################################################################
#
# Arc diagrams
#
################################################################################

const ArcDiagramVertex = Tuple{Symbol, Int}

abstract type ArcDiagram end

struct ArcDiagramUndirected <: ArcDiagram
    n_upper_verts::Int
    n_lower_verts::Int
    upper_neighbors::Vector{ArcDiagramVertex}
    lower_neighbors::Vector{ArcDiagramVertex}

    function ArcDiagramUndirected(
        n_upper_verts::Int,
        n_lower_verts::Int,
        upper_neighbors::Vector{ArcDiagramVertex},
        lower_neighbors::Vector{ArcDiagramVertex};
        check::Bool=true,
    )
        if check
            @req length(upper_neighbors) == n_upper_verts "Upper vertices' neighbors list has wrong length."
            @req length(lower_neighbors) == n_lower_verts "Lower vertices' neighbors list has wrong length."
            for i in 1:n_upper_verts
                neigh = upper_neighbors[i]
                if is_upper_vertex(neigh)
                    @req 1 <= vertex_index(neigh) <= n_upper_verts "Out of bounds adjacency."
                    @req vertex_index(neigh) != i "No self-loops allowed."
                    neigh_neigh = upper_neighbors[vertex_index(neigh)]
                    @req is_upper_vertex(neigh_neigh) && vertex_index(neigh_neigh) == i "Adjacency not symmetric."
                elseif is_lower_vertex(neigh)
                    @req 1 <= vertex_index(neigh) <= n_lower_verts "Out of bounds adjacency."
                    neigh_neigh = lower_neighbors[vertex_index(neigh)]
                    @req is_upper_vertex(neigh_neigh) && vertex_index(neigh_neigh) == i "Adjacency not symmetric."
                else
                    @req false "Invalid neighbor type."
                end
            end
            for i in 1:n_lower_verts
                neigh = lower_neighbors[i]
                if is_upper_vertex(neigh)
                    @req 1 <= vertex_index(neigh) <= n_upper_verts "Out of bounds adjacency."
                    neigh_neigh = upper_neighbors[vertex_index(neigh)]
                    @req is_lower_vertex(neigh_neigh) && vertex_index(neigh_neigh) == i "Adjacency not symmetric."
                elseif is_lower_vertex(neigh)
                    @req 1 <= vertex_index(neigh) <= n_lower_verts "Out of bounds adjacency."
                    @req vertex_index(neigh) != i "No self-loops allowed."
                    neigh_neigh = lower_neighbors[vertex_index(neigh)]
                    @req is_lower_vertex(neigh_neigh) && vertex_index(neigh_neigh) == i "Adjacency not symmetric."
                else
                    @req false "Invalid neighbor type."
                end
            end
        end
        return new(n_upper_verts, n_lower_verts, upper_neighbors, lower_neighbors)
    end

    function ArcDiagramUndirected(
        n_upper_verts::Int,
        n_lower_verts::Int,
        upper_neighbor_inds::Vector{Int},
        lower_neighbor_inds::Vector{Int};
        check::Bool=true,
    )
        upper_neighbors = [
            i < 0 ? upper_vertex(ArcDiagramVertex, -i) : lower_vertex(ArcDiagramVertex, i) for i in upper_neighbor_inds
        ]
        lower_neighbors = [
            i < 0 ? upper_vertex(ArcDiagramVertex, -i) : lower_vertex(ArcDiagramVertex, i) for i in lower_neighbor_inds
        ]
        return ArcDiagramUndirected(n_upper_verts, n_lower_verts, upper_neighbors, lower_neighbors; check)
    end
end

struct ArcDiagramDirected <: ArcDiagram
    n_upper_verts::Int
    n_lower_verts::Int
    parity_upper_verts::BitVector   # true is down
    parity_lower_verts::BitVector   # true is down
    upper_neighbors::Vector{ArcDiagramVertex}
    lower_neighbors::Vector{ArcDiagramVertex}

    function ArcDiagramDirected(
        n_upper_verts::Int,
        n_lower_verts::Int,
        parity_upper_verts::BitVector,
        parity_lower_verts::BitVector,
        upper_neighbors::Vector{ArcDiagramVertex},
        lower_neighbors::Vector{ArcDiagramVertex};
        check::Bool=true,
    )
        if check
            @req length(parity_upper_verts) == n_upper_verts "Upper vertices' parity list has wrong length."
            @req length(parity_lower_verts) == n_lower_verts "Lower vertices' parity list has wrong length."
            @req length(upper_neighbors) == n_upper_verts "Upper vertices' neighbors list has wrong length."
            @req length(lower_neighbors) == n_lower_verts "Lower vertices' neighbors list has wrong length."
            for i in 1:n_upper_verts
                neigh = upper_neighbors[i]
                if is_upper_vertex(neigh)
                    @req 1 <= vertex_index(neigh) <= n_upper_verts "Out of bounds adjacency."
                    @req vertex_index(neigh) != i "No self-loops allowed."
                    neigh_neigh = upper_neighbors[vertex_index(neigh)]
                    @req is_upper_vertex(neigh_neigh) && vertex_index(neigh_neigh) == i "Adjacency not symmetric."
                    @req parity_upper_verts[i] != parity_upper_verts[vertex_index(neigh)] "Direction parity mismatch."
                elseif is_lower_vertex(neigh)
                    @req 1 <= vertex_index(neigh) <= n_lower_verts "Out of bounds adjacency."
                    neigh_neigh = lower_neighbors[vertex_index(neigh)]
                    @req is_upper_vertex(neigh_neigh) && vertex_index(neigh_neigh) == i "Adjacency not symmetric."
                    @req parity_upper_verts[i] == parity_lower_verts[vertex_index(neigh)] "Direction parity mismatch."
                else
                    @req false "Invalid neighbor type."
                end
            end
            for i in 1:n_lower_verts
                neigh = lower_neighbors[i]
                if is_upper_vertex(neigh)
                    @req 1 <= vertex_index(neigh) <= n_upper_verts "Out of bounds adjacency."
                    neigh_neigh = upper_neighbors[vertex_index(neigh)]
                    @req is_lower_vertex(neigh_neigh) && vertex_index(neigh_neigh) == i "Adjacency not symmetric."
                    @req parity_lower_verts[i] == parity_upper_verts[vertex_index(neigh)] "Direction parity mismatch."
                elseif is_lower_vertex(neigh)
                    @req 1 <= vertex_index(neigh) <= n_lower_verts "Out of bounds adjacency."
                    @req vertex_index(neigh) != i "No self-loops allowed."
                    neigh_neigh = lower_neighbors[vertex_index(neigh)]
                    @req is_lower_vertex(neigh_neigh) && vertex_index(neigh_neigh) == i "Adjacency not symmetric."
                    @req parity_lower_verts[i] != parity_lower_verts[vertex_index(neigh)] "Direction parity mismatch."
                else
                    @req false "Invalid neighbor type."
                end
            end
        end
        return new(
            n_upper_verts,
            n_lower_verts,
            parity_upper_verts,
            parity_lower_verts,
            upper_neighbors,
            lower_neighbors,
        )
    end

    function ArcDiagramDirected(
        n_upper_verts::Int,
        n_lower_verts::Int,
        parity_upper_verts::BitVector,
        parity_lower_verts::BitVector,
        upper_neighbor_inds::Vector{Int},
        lower_neighbor_inds::Vector{Int};
        check::Bool=true,
    )
        upper_neighbors = [
            i < 0 ? upper_vertex(ArcDiagramVertex, -i) : lower_vertex(ArcDiagramVertex, i) for i in upper_neighbor_inds
        ]
        lower_neighbors = [
            i < 0 ? upper_vertex(ArcDiagramVertex, -i) : lower_vertex(ArcDiagramVertex, i) for i in lower_neighbor_inds
        ]
        return ArcDiagramDirected(
            n_upper_verts,
            n_lower_verts,
            parity_upper_verts,
            parity_lower_verts,
            upper_neighbors,
            lower_neighbors;
            check,
        )
    end
end

struct ArcDiagramIterator{T <: Union{Directed, Undirected}}
    iter
    len::Int
end

################################################################################
#
# Pseudographs
#
################################################################################

struct PseudographLabelled{T}
    nv::Int
    edges::MSet{Pair{MSet{Int}, T}}

    function PseudographLabelled(
        nv::Int,
        edges::MSet{Pair{MSet{Int}, T}};
        check::Bool=true,
        regular_degree::Union{Nothing, Int}=nothing,
    ) where {T}
        if check
            @req all(e -> all(i -> 1 <= i <= nv, first(e)), edges) "Out of bounds edge"
            if !isnothing(regular_degree)
                for i in 1:nv
                    @req regular_degree == sum(e -> multiplicity(first(e), i), edges; init=0) "Vertex $i has wrong degree"
                end
            end
        end
        return new{T}(nv, edges)
    end

    function PseudographLabelled(
        nv::Int,
        edges::Vector{Pair{MSet{Int}, T}};
        check::Bool=true,
        regular_degree::Union{Nothing, Int}=nothing,
    ) where {T}
        return PseudographLabelled(nv, MSet(edges); check, regular_degree)
    end
end

################################################################################
#
# Deformation maps
#
################################################################################

"""
    DeformationMap{C} = MatElem{FreeAssAlgElem{C}} where {C <: RingElem}

The type for deformation maps of a Lie algebra smash product.
The entry `kappa[i,j]` should be the image of ``v_i \\wedge v_j`` under the deformation map, i.e. ``κ(v_i,v_j)``.
Deformation maps are always assumed to be quadratic and skew-symmetric.
"""
const DeformationMap{C} = MatElem{<:FreeAssAlgElem{C}} where {C <: RingElem} # TODO: make concrete type


"""
    abstract type DeformBasis{C <: RingElem} end

A basis for a deformation map space of a Lie algebra smash product.
The constructor of a subtype should accept a [`SmashProductLie`](@ref) and an `AbstractVector{Int}` of degrees.
It is required that `Base.length` and `Base.iterate` are implemented for subtypes,
where iterating yields objects of type `DeformationMap{C}`.

For a reference implementation, we refer to [`StdDeformBasis`](@ref).
"""
abstract type DeformBasis{C <: RingElem} end

################################################################################
#
# Smash products
#
################################################################################

"""
The struct representing a Lie algebra smash product.
It consists of the underlying FreeAssAlgebra with relations and some metadata.
It gets created by calling [`smash_product`](@ref).
"""
@attributes mutable struct SmashProductLie{C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}} <: NCRing
    coeff_ring::Ring
    L::LieAlgebra{LieC}
    V::LieAlgebraModule{LieC}
    alg::FreeAssAlgebra{C}
    rels::Matrix{Union{Nothing, FreeAssAlgElem{C}}}

    # default constructor for @attributes
    function SmashProductLie{C, LieC, LieT}(
        coeff_ring::Ring,
        L::LieAlgebra{LieC},
        V::LieAlgebraModule{LieC},
        alg::FreeAssAlgebra{C},
        rels::Matrix{Union{Nothing, FreeAssAlgElem{C}}},
    ) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}}
        new{C, LieC, LieT}(coeff_ring, L, V, alg, rels)
    end
end

mutable struct SmashProductLieElem{C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}} <: NCRingElem
    p::SmashProductLie{C, LieC, LieT}   # parent
    alg_elem::FreeAssAlgElem{C}
    simplified::Bool

    function SmashProductLieElem(
        p::SmashProductLie{C, LieC, LieT},
        alg_elem::FreeAssAlgElem{C};
        simplified::Bool=false,
    ) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}}
        @req underlying_algebra(p) === parent(alg_elem) "Incompatible algebras."
        return new{C, LieC, LieT}(p, alg_elem, simplified)
    end
end

################################################################################
#
# Smash product deformations
#
################################################################################

"""
The struct representing a deformation of a Lie algebra smash product.
It consists of the underlying FreeAssAlgebra with relations and some metadata.
It gets created by calling [`deform`](@ref).
"""
@attributes mutable struct SmashProductLieDeform{C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}} <:
                           NCRing
    sp::SmashProductLie{C, LieC, LieT}
    rels::Matrix{Union{Nothing, FreeAssAlgElem{C}}}
    kappa::DeformationMap{C}

    # default constructor for @attributes
    function SmashProductLieDeform{C, LieC, LieT}(
        sp::SmashProductLie{C, LieC, LieT},
        rels::Matrix{Union{Nothing, FreeAssAlgElem{C}}},
        kappa::DeformationMap{C},
    ) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}}
        new{C, LieC, LieT}(sp, rels, kappa)
    end
end

mutable struct SmashProductLieDeformElem{C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}} <: NCRingElem
    p::SmashProductLieDeform{C, LieC, LieT}   # parent
    alg_elem::FreeAssAlgElem{C}
    simplified::Bool

    function SmashProductLieDeformElem(
        p::SmashProductLieDeform{C, LieC, LieT},
        alg_elem::FreeAssAlgElem{C};
        simplified::Bool=false,
    ) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}}
        @req underlying_algebra(p) === parent(alg_elem) "Incompatible algebras."
        return new{C, LieC, LieT}(p, alg_elem, simplified)
    end
end
