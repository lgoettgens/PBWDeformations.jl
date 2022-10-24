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
Given a basis element of an above basis, there are functions to find one/all arc diagrams that induce it (up to a scalar). This is currently done by recomputing the basis and keeping track of the used arc diagrams.
```@docs
corresponding_arc_diagram
corresponding_arc_diagrams
```
