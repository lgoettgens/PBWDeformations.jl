using Combinatorics
using CSV
using DataFrames

function summarizeDB(path, delim)
    inputs = unique(split(filename, '_')[1:3] for filename in readdir(path) if filename[2] == '_')

    rows = Tuple{Char, Int, Vector{Int}, Int, Vector{Int}}[]

    for input in inputs
        dynkin = input[1][1]

        n = parse(Int, input[2])

        regex = r"([0-9])"
        lambda = [parse(Int, t.match) for t in eachmatch(regex, input[3])]
        lambda_string = replace("$lambda", ' ' => "")

        prefix = "$(path)/$(dynkin)_$(n)_$(lambda_string)_"
        suffix = ".txt"
        degs = map(filename -> parse(Int, filename[length(prefix)+1:end-length(suffix)]), filter(startswith(prefix), readdir(path; join=true)))
        if isempty(degs)
            @warn "No valid file for ($(dynkin), $(n), $(lambda_string)) found."
            continue
        end
        maxdeg = maximum(degs)
        if !issetequal(degs, 0:maxdeg)
            @warn "Files for ($(dynkin), $(n), $(lambda_string)) are missing some degree(s)."
            continue
        end
        new_deform_degs = Int[]
        for deg in degs
            if deg == 0 && filesize("$(prefix)$(deg)$(suffix)") > 3
                push!(new_deform_degs, deg)
            end
            if deg > 0 && filesize("$(prefix)$(deg)$(suffix)") > filesize("$(prefix)$(deg-1)$(suffix)")
                push!(new_deform_degs, deg)
            end
        end
        push!(rows, (dynkin, n, lambda, maxdeg, new_deform_degs))
    end

    CSV.write(
        "$(path)/summary.csv",
        DataFrame(NamedTuple{(:dynkin, :n, :lambda, :maxdeg, :new_deform_degs)}.(rows));
        delim=delim,
        quotestrings=false,
        transform=(col, val) -> (typeof(val) <: Vector && isempty(val) ? "[]" : replace("$val", ' ' => "")),
    )
end

path = "db"
csvdelim = '\t'

summarizeDB(path, csvdelim)
