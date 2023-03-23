@info "Parsing ARGS..."

@assert length(ARGS[1]) == 1 || ARGS[1] == "SO"
if length(ARGS[1]) == 1
    dynkin = ARGS[1][1]
else
    dynkin = "SO"
end
n = parse(Int, ARGS[2])

regex = r"([0-9])"
lambda = [parse(Int, t.match) for t in eachmatch(regex, ARGS[3])]
lambda_string = replace("$lambda", ' ' => "")

maxdeg = parse(Int, ARGS[4])

if isa(dynkin, Char)
    @assert n == length(lambda)
elseif dynkin == "SO"
    @assert div(n, 2) == length(lambda)
else
    @assert false
end

@info "Loading Dependencies..."
using Oscar
using PBWDeformations

function kappa_to_string(kappa::Matrix{FreeAssAlgElem{fmpq}})
    dim = size(kappa, 1)
    string = "["
    for i in 1:dim
        for j in 1:dim
            string *= "("
            string *= sprint(show, kappa[i, j])
            string *= ")"
            if j < dim
                string *= "\t"
            end
        end
        if i == dim
            string *= "]"
        end
        string *= "\n"
    end
    return string
end

for deg in 0:maxdeg
    @info "Computing PBWDeformations... (degree=$deg)"
    if isa(dynkin, Char)
        sp, _ = smash_product_lie(QQ, dynkin, n, lambda)
    elseif dynkin == "SO"
        sp, _ = smash_product_lie_so(QQ, n, lambda)
    else
        @assert false
    end
    kappas = all_pbwdeformations(sp, deg)

    @info "Writing Output... (degree=$deg)"
    path = "db"
    filename = "$(path)/$(dynkin)_$(n)_$(lambda_string)_$(deg).txt"
    Base.Filesystem.mkpath(path)
    open(filename, "w") do io
        for kappa in kappas
            println(io, kappa_to_string(kappa))
        end
        print(io, "EOF")
    end
end
