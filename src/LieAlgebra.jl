@attributes mutable struct LieAlgebra{C <: RingElement} <: FPModule{C}
    R::Ring
    n::Int
    basis::Vector{<:MatElem{C}}
    dim::Int
    s::Vector{Symbol}

    function LieAlgebra{C}(
        R::Ring,
        n::Int,
        basis::Vector{<:MatElem{C}},
        s::Vector{Symbol},
        cached::Bool=true,
    ) where {C <: RingElement}
        return get_cached!(LieAlgebraDict, (R, n, basis, s), cached) do
            all(b -> size(b) == (n, n), basis) || error("Invalid basis element dimensions.")
            length(s) == length(basis) || error("Invalid number of basis element names.")
            new{C}(R, n, basis, length(basis), s)
        end::LieAlgebra{C}
    end

    function LieAlgebra{C}(
        R::Ring,
        n::Int,
        basis::Vector{<:MatElem{C}},
        s::Vector{<:Union{AbstractString, Char, Symbol}},
        cached::Bool=true,
    ) where {C <: RingElement}
        return LieAlgebra{C}(R, n, basis, Symbol.(s), cached)
    end
end

const LieAlgebraDict = CacheDictType{Tuple{Ring, Int, Vector{<:MatElem}, Vector{Symbol}}, LieAlgebra}()

struct LieAlgebraElem{C <: RingElement} <: FPModuleElem{C}
    parent::LieAlgebra{C}
    mat::MatElem{C}
end


###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{LieAlgebraElem{C}}) where {C <: RingElement} = LieAlgebra{C}

elem_type(::Type{LieAlgebra{C}}) where {C <: RingElement} = LieAlgebraElem{C}

parent(x::LieAlgebraElem{C}) where {C <: RingElement} = x.parent

base_ring(L::LieAlgebra{C}) where {C <: RingElement} = L.R

base_ring(x::LieAlgebraElem{C}) where {C <: RingElement} = base_ring(parent(x))

ngens(L::LieAlgebra{C}) where {C <: RingElement} = L.dim

gens(L::LieAlgebra{C}) where {C <: RingElement} = [gen(L, i) for i in 1:ngens(L)]

function gen(L::LieAlgebra{C}, i::Int) where {C <: RingElement}
    R = base_ring(L)
    return L([(j == i ? one(R) : zero(R)) for j in 1:ngens(L)])
end

@inline function Generic._matrix(x::LieAlgebraElem{T}) where {T}
    return (x.mat)::dense_matrix_type(T)
end

@inline function Generic.basis(L::LieAlgebra{T}) where {T}
    return (L.basis)::Vector{dense_matrix_type(T)}
end

function Generic.rels(_::LieAlgebra{C}) where {C <: RingElement}
    # there are no relations in a vector space
    return Vector{dense_matrix_type(C)}(undef, 0)
end

###############################################################################
#
#   String I/O
#
###############################################################################

# TODO

function Base.show(io::IO, x::LieAlgebraElem{C}) where {C <: RingElement}
    Base.show(io, matrix_repr(x))
end


###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (L::LieAlgebra{C})() where {C <: RingElement}
    mat = zero_matrix(base_ring(L), 1, ngens(L))
    return elem_type(L)(L, mat)
end

function (L::LieAlgebra{C})(v::Vector{C}) where {C <: RingElement}
    length(v) == ngens(L) || error("Length of vector does not match number of generators.")
    mat = matrix(base_ring(L), 1, length(v), v)
    return elem_type(L)(L, mat)
end

function (L::LieAlgebra{C})(m::MatElem{C}) where {C <: RingElement}
    if size(m) == (L.n, L.n)
        m = coefficient_vector(m, basis(L))
    end
    size(m) == (1, ngens(L)) || error("Invalid matrix dimensions.")
    return elem_type(L)(L, m)
end

function (L::LieAlgebra{C})(v::LieAlgebraElem{C}) where {C <: RingElement}
    L == parent(v) || error("Incompatible modules.")
    return v
end


###############################################################################
#
#   Arithmetic operations
#
###############################################################################

# Vector space operations get inherited from FPModule

function Generic.matrix_repr(x::LieAlgebraElem{C}) where {C <: RingElement}
    return sum(c * b for (c, b) in zip(x.mat, basis(parent(x))))
end

function bracket(x::LieAlgebraElem{C}, y::LieAlgebraElem{C}) where {C <: RingElement}
    check_parent(x, y)
    L = parent(x)
    x_mat = matrix_repr(x)
    y_mat = matrix_repr(y)
    return L(x_mat * y_mat - y_mat * x_mat)
end


###############################################################################
#
#   Comparison functions
#
###############################################################################

# Overwrite the equality of FPModule to be used for CacheDicts
function Base.:(==)(L1::LieAlgebra{C}, L2::LieAlgebra{C}) where {C <: RingElement}
    return L1 === L2
end


###############################################################################
#
#   Constructor
#
###############################################################################

function general_linear_liealgebra(R::Ring, n::Int)
    basis = [(b = zero_matrix(R, n, n); b[i, j] = 1; b) for i in 1:n for j in 1:n]
    s = ["x_$(i)_$(j)" for i in 1:n for j in 1:n]
    L = LieAlgebra{elem_type(R)}(R, n, basis, s)
    set_attribute!(L, :type, :general_linear)
    return L
end

function special_linear_liealgebra(R::Ring, n::Int)
    basis_e = [(b = zero_matrix(R, n, n); b[i, j] = 1; b) for i in 1:n for j in i+1:n]
    basis_f = [(b = zero_matrix(R, n, n); b[j, i] = -1; b) for i in 1:n for j in i+1:n]
    basis_h = [(b = zero_matrix(R, n, n); b[i, i] = 1; b[i+1, i+1] = -1; b) for i in 1:n-1]
    s_e = ["e_$(i)_$(j)" for i in 1:n for j in i+1:n]
    s_f = ["f_$(i)_$(j)" for i in 1:n for j in i+1:n]
    s_h = ["h_$(i)" for i in 1:n-1]
    L = LieAlgebra{elem_type(R)}(R, n, [basis_e; basis_f; basis_h], [s_e; s_f; s_h])
    set_attribute!(L, :type, :special_linear)
    return L
end

function special_orthogonal_liealgebra(R::Ring, n::Int)
    basis = [(b = zero_matrix(R, n, n); b[i, j] = 1; b[j, i] = -1; b) for i in 1:n for j in i+1:n]
    s = ["x_$(i)_$(j)" for i in 1:n for j in i+1:n]
    L = LieAlgebra{elem_type(R)}(R, n, basis, s)
    set_attribute!(L, :type, :special_orthogonal)
    return L
end
