mutable struct LieAlgebraAbstractModule{C <: RingElement} <: LieAlgebraModule{C}
    L::LieAlgebra{C}
    dim::Int
    transformation_matrices::Vector{MatElem{C}}
    s::Vector{Symbol}

    function LieAlgebraAbstractModule{C}(
        L::LieAlgebra{C},
        dim::Int,
        transformation_matrices::Vector{<:MatElem{C}},
        s::Vector{Symbol},
        cached::Bool=true;
        check::Bool=true,
    ) where {C <: RingElement}
        return get_cached!(
            LieAlgebraAbstractModuleDict,
            (L, dim, transformation_matrices, s),
            cached,
        ) do
            all(m -> size(m) == (dim, dim), transformation_matrices) ||
                error("Invalid transformation matrix dimensions.")
            dim == length(s) || error("Invalid number of basis element names.")

            V = new{C}(L, dim, transformation_matrices, s)
            if check
                for xi in gens(L), xj in gens(L), v in gens(V)
                    bracket(xi, xj) * v == xi * (xj * v) - xj * (xi * v) ||
                        error("Structure constants do not define a module.")
                end
            end
            V
        end::LieAlgebraAbstractModule{C}
    end

end

const LieAlgebraAbstractModuleDict =
    CacheDictType{Tuple{LieAlgebra, Int, Vector{MatElem}, Vector{Symbol}}, LieAlgebraAbstractModule}()

struct LieAlgebraAbstractModuleElem{C <: RingElement} <: LieAlgebraModuleElem{C}
    parent::LieAlgebraAbstractModule{C}
    mat::MatElem{C}
end


###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{LieAlgebraAbstractModuleElem{C}}) where {C <: RingElement} = LieAlgebraAbstractModule{C}

elem_type(::Type{LieAlgebraAbstractModule{C}}) where {C <: RingElement} = LieAlgebraAbstractModuleElem{C}

parent(v::LieAlgebraAbstractModuleElem{C}) where {C <: RingElement} = v.parent

base_ring(V::LieAlgebraAbstractModule{C}) where {C <: RingElement} = base_ring(base_liealgebra(V))

base_liealgebra(V::LieAlgebraAbstractModule{C}) where {C <: RingElement} = V.L

ngens(V::LieAlgebraAbstractModule{C}) where {C <: RingElement} = V.dim


###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, V::LieAlgebraAbstractModule{C}) where {C <: RingElement}
    print(io, "AbstractModule of ")
    print(IOContext(io, :compact => true), base_liealgebra(V))
end

function symbols(V::LieAlgebraAbstractModule{C}) where {C <: RingElement}
    return V.s
end


###############################################################################
#
#   Parent object call overload
#
###############################################################################

# no special ones


###############################################################################
#
#   Module action
#
###############################################################################

function transformation_matrix_by_basisindex(V::LieAlgebraAbstractModule{C}, i::Int) where {C <: RingElement}
    return (V.transformation_matrices[i])::dense_matrix_type(C)
end


###############################################################################
#
#   Constructor
#
###############################################################################

function abstract_module(
    L::LieAlgebra{C},
    dim::Int,
    transformation_matrices::Vector{<:MatElem{C}},
    s::Vector{<:Union{AbstractString, Char, Symbol}},
    cached::Bool=true;
    check::Bool=true,
) where {C <: RingElement}
    return LieAlgebraAbstractModule{C}(L, dim, transformation_matrices, Symbol.(s), cached; check=check)
end

function abstract_module(
    L::LieAlgebra{C},
    dim::Int,
    struct_consts::Matrix{SRow{C}},
    s::Vector{<:Union{AbstractString, Char, Symbol}},
    cached::Bool=true;
    check::Bool=true,
) where {C <: RingElement}
    ngens(L) == size(struct_consts, 1) || error("Invalid structure constants dimensions.")
    dim == size(struct_consts, 2) || error("Invalid structure constants dimensions.")
    dim == length(s) || error("Invalid number of basis element names.")

    transformation_matrices = [zero_matrix(base_ring(L), dim, dim) for _ in 1:ngens(L)]
    for i in 1:ngens(L), j in 1:dim
        transformation_matrices[i][:, j] = transpose(dense_row(struct_consts[i, j], dim))
    end

    return LieAlgebraAbstractModule{C}(L, dim, transformation_matrices, Symbol.(s), cached; check=check)
end
