using Oscar
import AbstractAlgebra: RingElement

module DenseMatrixColMajor

using SparseArrays
using Oscar

import ..RingElement

function findfirstnz(v::Vector{T}, start=1::Int) where {T <: Union{RingElement, Number}}
    for i in start:length(v)
        if !iszero(v[i])
            return i, v[i]
        end
    end
    return 0, 0
end

function reduce_and_store!(lgs::Matrix{T}, v::Vector{T}) where {T <: Union{RingElement, Number}}
    nz_ind = 0
    while !iszero(v)
        nz_ind, nz_val = findfirstnz(v, nz_ind + 1)
        v = inv(nz_val) .* v
        if iszero(lgs[nz_ind, nz_ind])
            lgs[:, nz_ind] = v
            return
        else
            v -= lgs[:, nz_ind]
        end
    end
end

function reduced_row_echelon!(lgs::Matrix{T}) where {T <: Union{RingElement, Number}}
    for i in size(lgs, 2):-1:1
        if iszero(lgs[i, i])
            continue
        end
        nz_inds, nz_vals = findnz(sparse(lgs[:, i]))
        for (ind, j) in enumerate(nz_inds[2:end])
            if !iszero(lgs[j, j])
                lgs[:, i] -= nz_vals[ind+1] .* lgs[:, j]
            end
        end
    end
    return lgs
end

function lgs_to_mat(lgs::Matrix{T}) where {T <: Union{RingElement, Number}}
    return sparse(lgs)
end

function indices_of_freedom(mat::SparseArrays.SparseMatrixCSC{T, Int}) where {T <: Union{RingElement, Number}}
    size(mat)[1] == size(mat)[2] || throw(ArgumentError("Matrix needs to be square."))
    return filter(i -> iszero(mat[i, i]), 1:size(mat)[1])
end

function solve_lgs(instance)
    iter = instance

    lgs = zeros(typeof(instance[1][1]), length(instance[1]), length(instance[1]))

    for v in iter
        reduce_and_store!(lgs, v)
    end

    reduced_row_echelon!(lgs)

    mat = lgs_to_mat(lgs)

    freedom_ind = indices_of_freedom(mat)
    freedom_deg = length(freedom_ind)

    return freedom_deg, freedom_ind
end

end

module Dense

using SparseArrays
using Oscar

import ..RingElement

function findfirstnz(v::Vector{T}, start=1::Int) where {T <: Union{RingElement, Number}}
    for i in start:length(v)
        if !iszero(v[i])
            return i, v[i]
        end
    end
    return 0, 0
end

function reduce_and_store!(lgs::Vector{Union{Nothing, Vector{T}}}, v::Vector{T}) where {T <: Union{RingElement, Number}}
    nz_ind = 0
    while !iszero(v)
        nz_ind, nz_val = findfirstnz(v, nz_ind + 1)
        v = inv(nz_val) .* v
        if lgs[nz_ind] === nothing
            lgs[nz_ind] = v
            return
        else
            v -= lgs[nz_ind]
        end
    end
end

function reduced_row_echelon!(lgs::Vector{Union{Nothing, Vector{T}}}) where {T <: Union{RingElement, Number}}
    for i in length(lgs):-1:1
        if lgs[i] === nothing
            continue
        end
        nz_inds, nz_vals = findnz(sparse(lgs[i]))
        for (ind, j) in enumerate(nz_inds[2:end])
            if lgs[j] !== nothing
                lgs[i] -= nz_vals[ind+1] .* lgs[j]
            end
        end
    end
    return lgs
end

function lgs_to_mat(lgs::Vector{Union{Nothing, Vector{T}}}) where {T <: Union{RingElement, Number}}
    n = length(lgs)
    mat = spzeros(T, n, n)
    for i in 1:n
        if lgs[i] !== nothing
            mat[i, :] = lgs[i]
        end
    end
    return mat
end

function indices_of_freedom(mat::SparseArrays.SparseMatrixCSC{T, Int}) where {T <: Union{RingElement, Number}}
    size(mat)[1] == size(mat)[2] || throw(ArgumentError("Matrix needs to be square."))
    return filter(i -> iszero(mat[i, i]), 1:size(mat)[1])
end

function solve_lgs(instance)
    iter = instance

    lgs = Vector{Union{Nothing, Vector{fmpq}}}(nothing, length(instance[1]))
    for v in iter
        reduce_and_store!(lgs, v)
    end

    reduced_row_echelon!(lgs)

    mat = lgs_to_mat(lgs)

    freedom_ind = indices_of_freedom(mat)
    freedom_deg = length(freedom_ind)

    return freedom_deg, freedom_ind
end

end

module SparseCurrent

using SparseArrays
using Oscar

import ..RingElement

function reduce_and_store!(
    lgs::Vector{Union{Nothing, SparseVector{T, Int}}},
    v::SparseVector{T, Int},
) where {T <: Union{RingElement, Number}}
    while !iszero(v)
        nz_inds, nz_vals = findnz(v)
        v = inv(nz_vals[1]) .* v
        if lgs[nz_inds[1]] === nothing
            lgs[nz_inds[1]] = v
            return
        else
            v -= lgs[nz_inds[1]]
        end
    end
end

function reduced_row_echelon!(lgs::Vector{Union{Nothing, SparseVector{T, Int}}}) where {T <: Union{RingElement, Number}}
    for i in length(lgs):-1:1
        if lgs[i] === nothing
            continue
        end
        nz_inds, nz_vals = findnz(lgs[i])
        for (ind, j) in enumerate(nz_inds[2:end])
            if lgs[j] !== nothing
                lgs[i] -= nz_vals[ind+1] .* lgs[j]
            end
        end
    end
    return lgs
end

function lgs_to_mat(lgs::Vector{Union{Nothing, SparseVector{T, Int}}}) where {T <: Union{RingElement, Number}}
    n = length(lgs)
    mat = spzeros(T, n, n)
    for i in 1:n
        if lgs[i] !== nothing
            mat[i, :] = lgs[i]
        end
    end
    return mat
end

function indices_of_freedom(mat::SparseArrays.SparseMatrixCSC{T, Int}) where {T <: Union{RingElement, Number}}
    size(mat)[1] == size(mat)[2] || throw(ArgumentError("Matrix needs to be square."))
    return filter(i -> iszero(mat[i, i]), 1:size(mat)[1])
end

function solve_lgs(instance)
    iter = Iterators.map(sparse, instance)

    lgs = Vector{Union{Nothing, SparseVector{fmpq, Int}}}(nothing, length(instance[1]))
    for v in iter
        reduce_and_store!(lgs, v)
    end

    reduced_row_echelon!(lgs)

    mat = lgs_to_mat(lgs)

    freedom_ind = indices_of_freedom(mat)
    freedom_deg = length(freedom_ind)

    return freedom_deg, freedom_ind
end

end

module SparseBetter

using SparseArrays
using Oscar

import ..RingElement

function reduce_and_store!(
    lgs::Vector{Union{Nothing, SparseVector{T, Int}}},
    v::SparseVector{T, Int},
) where {T <: Union{RingElement, Number}}
    while count(!iszero, v) > 0
        nz_inds, nz_vals = findnz(v)
        if !isone(nz_vals[1])
            v = inv(nz_vals[1]) .* v
        end
        if lgs[nz_inds[1]] === nothing
            lgs[nz_inds[1]] = v
            return
        else
            v -= lgs[nz_inds[1]]
        end
    end
end

function reduced_row_echelon!(lgs::Vector{Union{Nothing, SparseVector{T, Int}}}) where {T <: Union{RingElement, Number}}
    for i in length(lgs):-1:1
        if lgs[i] === nothing
            continue
        end
        nz_inds, nz_vals = findnz(lgs[i])
        for (ind, j) in enumerate(nz_inds[2:end])
            if lgs[j] !== nothing
                lgs[i] -= nz_vals[ind+1] .* lgs[j]
            end
        end
    end
    return lgs
end

function lgs_to_mat(lgs::Vector{Union{Nothing, SparseVector{T, Int}}}) where {T <: Union{RingElement, Number}}
    n = length(lgs)
    mat = spzeros(T, n, n)
    for i in 1:n
        if lgs[i] !== nothing
            mat[i, :] = lgs[i]
        end
    end
    return mat
end

function indices_of_freedom(mat::SparseArrays.SparseMatrixCSC{T, Int}) where {T <: Union{RingElement, Number}}
    size(mat)[1] == size(mat)[2] || throw(ArgumentError("Matrix needs to be square."))
    return filter(i -> iszero(mat[i, i]), 1:size(mat)[1])
end

function solve_lgs(instance)
    iter = Iterators.map(sparse, instance)

    lgs = Vector{Union{Nothing, SparseVector{fmpq, Int}}}(nothing, length(instance[1]))
    for v in iter
        reduce_and_store!(lgs, v)
    end

    reduced_row_echelon!(lgs)

    mat = lgs_to_mat(lgs)

    freedom_ind = indices_of_freedom(mat)
    freedom_deg = length(freedom_ind)

    return freedom_deg, freedom_ind
end

end
