```@meta
CurrentModule = PBWDeformations
DocTestSetup  = quote
    using PBWDeformations
    using Oscar
end
```

# Pseudographs

```@autodocs
Modules = [PBWDeformations]
Pages   = ["Pseudograph.jl"]
```

## Pseudograph induced bases
!!! warning
    The basis [`PseudographDeformBasis`](@ref) can currently only be used for exterior powers of the standard module of Lie type $\mathfrak{so}_n$.
```@docs
PseudographDeformBasis
```

## Reverse direction
Given a basis element of an above basis, one can lookup all pseudographs that induce it (up to a scalar).
See [`lookup_data`](@ref) for more details.
