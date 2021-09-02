
sol = nullspace(Matrix(A)) # here, we solve the LSE
xmax(coll) = isempty(coll) ? 0 : maximum(coll)
if xmax(abs.(A*sol)) > 1e-6 warn("numerical imprecision") end

# reformat solutions...
sol = sparse(sol)
#display(sol)
say("$(sol.n) solutions:")
flush(stdout)

solEx = []
for j = 1:sol.n
v = sol[:,j]
av = abs.(v)
v ./= minimum(av[av.>=1e-6])
v = (x->rationalize(x, tol=1e-4)).(v)
if sum(v.<0) > sum(v.>0) v=-v end
#v = myrat(v, factorial(maximum(dgs))*2^(maximum(dgs)))
#v = (x->round(x, sigdigits=6)).(v)
push!(solEx, sum( v[i]*x(m...) for (i,m) in enumerate(mons) ))
end
