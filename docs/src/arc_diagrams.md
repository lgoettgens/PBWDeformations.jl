```@meta
CurrentModule = PBWDeformations
DocTestSetup  = quote
    using PBWDeformations
    using Oscar
end
```

# Arc diagrams

```@autodocs
Modules = [PBWDeformations]
Pages   = ["ArcDiagram.jl"]
```

## Arc diagram induced bases
!!! warning
    The basis [`ArcDiagDeformBasis`](@ref) can currently only be used for exterior powers of the standard module of Lie type $\mathfrak{so}_n$.
```@docs
ArcDiagDeformBasis
```

## Reverse direction
Given a basis element of an above basis, one can lookup all arc diagrams that induce it (up to a scalar).
```@docs
lookup_data
```
