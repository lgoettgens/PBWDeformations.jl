using Combinatorics
using Distributed

timeout = "60s"

max_n = 5
max_lambdasum = 3
maxdeg = 3
inputs = Tuple{Union{Char, String}, Int64, Vector{Int64}, Int64}[]

function lambdas(n, lambdasum)
    multisets(n, k) = map(A -> [sum(A .== i) for i in 1:n], with_replacement_combinations(1:n, k))
    return multisets(n, lambdasum)
end

#! format: off

# dynkin A
append!(inputs, [('A', n, lambda, maxdeg)
    for n in 1:max_n
    for lambdasum in 1:max_lambdasum
    for lambda in lambdas(n, lambdasum)])

# dynkin B and C
append!(inputs, [('B', n, lambda, maxdeg)
    for n in 2:max_n
    for lambdasum in 1:max_lambdasum
    for lambda in lambdas(n, lambdasum)])
append!(inputs, [('C', n, lambda, maxdeg)
    for n in 2:max_n
    for lambdasum in 1:max_lambdasum
    for lambda in lambdas(n, lambdasum)])

# dynkin D
append!(inputs, [('D', n, lambda, maxdeg)
    for n in 4:max_n
    for lambdasum in 1:max_lambdasum
    for lambda in lambdas(n, lambdasum)])

# other SO version
append!(inputs, [("SO", n, lambda, maxdeg)
    for n in 1:max_n
    for lambdasum in 1:max_lambdasum
    for lambda in lambdas(div(n,2), lambdasum)])

#! format: on

@info "Computing $(length(inputs)) cases..."
runtimes = pmap(inputs) do input
    println(input)
    dynkin, n, lambda, maxdeg = input
    lambda_string = replace("$lambda", ' ' => "")
    command = Cmd(
        `timeout -s SIGKILL $(timeout) julia scripts/CreateDBEntry.jl $(dynkin) $(n) $(lambda_string) $(maxdeg)`,
    )
    command = Cmd(command, ignorestatus=true)
    runtime = @elapsed run(command)
    return round(Int, runtime)
end

@info "Writing runtimes to file..."
path = "db"
filename = "$(path)/_runtimes.txt"
Base.Filesystem.mkpath(path)
open(filename, "w") do io
    for (input, runtime) in zip(inputs, runtimes)
        println(io, "$(input) => $(runtime)s")
    end
    print(io, "EOF")
end
