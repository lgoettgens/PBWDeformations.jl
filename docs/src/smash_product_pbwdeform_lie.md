```@meta
CurrentModule = PBWDeformations
DocTestSetup  = quote
    using PBWDeformations
    using Oscar
end
```

# PBW Deformations of smash products

## General deformation functions

```@docs
is_pbwdeformation
pbwdeform_eqs
```

## All PBW deformations

```@docs
all_pbwdeformations
```

### Bases of deformation map spaces

```@docs
DeformBasis
```

#### Standard basis

```@docs
StdDeformBasis
```

#### Other bases

Please refer to [Arc diagram induced bases](@ref) and [Pseudograph induced bases](@ref) for more specialized bases.
