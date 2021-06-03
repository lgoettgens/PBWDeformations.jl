
module SmashSO

using SymPy, Combinatorics

n = -1
comm_table = Dict()
x = sympy.Function("x", commutative=false)
v(i) = x(i)

base() = n < 10 ? 10 : n  # use base 10 for small n for readible expressions
q(i::Int,j::Int)::Int = base()*i+j

function init(_n)
  global n = _n
  
  ct = Dict()
  for i = 1:n, j=1:n, k=1:n
      if k in (i,j)
        ct[(q(i,j), k)] = j==k ? [(1, i)] : [(-1, j)]
      end
  end
  
  for i = 1:n, j = i+1:n, k = 1:n, l = k+1:n
      if q(i,j) <= q(k,l) continue end
      res = Array{Tuple{Int,Int}}([])
      add(f,i,j) = push!(res, i<j ? (f, q(i,j)) : (-f, q(j,i)))
      if j == k add(1,i,l) end
      if j == l add(-1,i,k) end
      if i == k add(-1,j,l) end
      if i == l add(1,j,k) end
      if res != [] ct[(q(i,j), q(k,l))] = res end
  end
  #display(ct)
  global comm_table = ct
  return
end

init(3) # default

function _normalForm(ind, ff=1)
    i = findfirst(ind[1:end-1] .> ind[2:end])
    if i == nothing return [(ff, ind)] end
    p = (ind[i], ind[i+1])
    commutators = haskey(comm_table, p) ? comm_table[p] : []
    
    return [
      _normalForm([ind[1:i-1]..., ind[i+1], ind[i], ind[i+2:end]...], ff)...,
      vcat((_normalForm([ind[1:i-1]..., j, ind[i+2:end]...], ff*f)
        for (f,j) in commutators)...)...
    ]
end

function normalForm(ind)
  xsum(coll) = isempty(coll) ? 0 : sum(coll)
  return xsum(f*x(nind...) for (f,nind) in _normalForm(ind))
end


normal(ex) = ex.replace(f -> f.func == x,
  f -> normalForm([N.(f.args)...])
)
comm(ex, i::Int) = ex.replace(f -> f.func == x,
  f -> normalForm([N.(f.args)..., i])
) - ex.replace(f -> f.func == x,
  f -> normalForm([i, N.(f.args)...])
).expand()


function test()
  println("testing...")
  @assert n >= 3
  @assert comm(x(1), 2) == 0
  @assert comm(x(12), 2) == x(1)
  @assert comm(x(12), 1) == -x(2)
  @assert comm(x(12), 13) == -x(23)
  one = sympify(1)
  @assert comm(one, 1) == 0
  @assert comm(one, 12) == 0
  @assert normalForm([2, 1]) == x(1, 2)
  #display(normalForm([12, 12, 2]))
  @assert normalForm([12, 12, 2]) == x(2,12,12) + 2*x(1,12) - x(2)
  #display(comm(x(2, 12), 2))
  @assert comm(x(2, 12), 2) == x(1, 2)
  println("done.")
  flush(stdout)
end

test()


function makeEl(dv, dg)
  if n < 1 error("n < 1 ?!") end
  pairs = [ q(i,j) for i=1:n for j=i+1:n ]
  monomials = [ x(itv..., itg...)
      for itv in powerset(1:n, 0, dv)
      for itg in powerset(pairs, 0, dg)
  ]
  return sum( Sym("a_$i")*mon for (i,mon) in enumerate(monomials) )
end

#v = sympy.Sym("v")
letter(i, nc=true) = i < base() ? sympy.Symbol("v_$i", commutative=!nc) : sympy.Symbol("X_$i", commutative=!nc)
format(ex, nc=true) = ex.replace(f -> f.func == x,
  f -> prod(letter(i, nc) for i in N.(f.args))
)
vars() = [(Sym("v_$k") for k in 1:n)..., (Sym("X_$k") for k in [q(i,j) for i=1:n for j=i+1:n])...]
coeffs(ex) = sympy.poly(format(ex, false), vars()...).coeffs()
show(ex) = (display.(format.(ex)); return)

using SparseArrays, LinearAlgebra

warn(msg) = (println("!!! $msg"); flush(stdout))



#display(vars())



### simfy
function _collect(factors)
  factors = [factors...]
  n = length(factors)
  for i=1:n-1
    if factors[i].func == x && factors[i+1].func == x
      factors[i+1] = x(factors[i].args..., factors[i+1].args...)
      factors[i] = 1
    end
  end
  return normal(prod(factors))
end

simfy(ex) = ex.expand().replace(sympy.Pow,
  (el, p) -> _collect([el for _=1:N(p)])
).replace(sympy.Mul,
  (m...) -> _collect(m)
).expand()

function test_simfy()
  @assert simfy(x(1, 2)*x(3, 4)) == x(1, 2, 3, 4)
  @assert simfy(x(12)*x(1, 2)) == -x(2, 2) + x(1, 1) + x(1, 2, 12)
end

test_simfy()

#t2 = simfy(t*t)
#show(t2)

end # module
