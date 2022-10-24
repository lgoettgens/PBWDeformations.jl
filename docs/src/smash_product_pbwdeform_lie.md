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
is_pbwdeform
pbwdeform_eqs
```

## All PBW deformations

```@docs
pbwdeforms_all
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

Please refer to [Arc diagram induced bases](@ref) for more specialized bases.
