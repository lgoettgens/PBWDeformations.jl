using Pkg
rev = readline("OSCAR.rev")
Pkg.add(PackageSpec(; url="https://github.com/oscar-system/Oscar.jl", rev))
