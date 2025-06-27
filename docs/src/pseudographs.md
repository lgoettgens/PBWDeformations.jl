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
    The basis [`PseudographDeformBasis`](@ref) can currently only be used for exterior and symmetric powers of the standard module of special orthogonal Lie algebras.
```@docs
PseudographDeformBasis
```

## Reverse direction
Given a basis element of an above basis, one can lookup all pseudographs that induce it (up to a scalar).
See [`lookup_params`](@ref) for more details.
