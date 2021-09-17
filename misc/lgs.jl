function posInMap!(map::Dict, el)
    if !(el in keys(map)) map[el] = map.count+1 end
    return map[el]
end

function sameLength(vecs...)
     n = maximum(v.n for v in vecs)
     return ( sparsevec(v.nzind, v.nzval, n) for v in vecs )    
end


  map = Dict{Vector{Int},Int}()
  indToPos(ind::Vector{Int}) = posInMap!(map, ind)
  collToVec(c::Vector{Tuple{Int, Vector{Int}}}) =
    sparsevec(indToPos.(last.(c)), first.(c)) # default: combine with +
  
  
  # set up A, a matrix encoding a linear system of equations
  mats = []
  
  for i in nums
    vecs = SparseVector[]

    a = collToVec(...)
    push!(vecs, a)

    end
    
  A = hcat(sameLength(vecs...)...)

  dropzeros!(A)
  rows = unique(rowvals(A))
  A = A[rows,:]
  say("($(size(A)[1]) eqs, $(size(A)[2]) vars)")
  if verbose display(A) end
  
  sol = nullspace(Matrix(A)) # here, we solve the LSE
  xmax(coll) = isempty(coll) ? 0 : maximum(coll)
  if xmax(abs.(A*sol)) > 1e-6 warn("numerical imprecision") end

  # reformat solutions...
  sol = sparse(sol)


