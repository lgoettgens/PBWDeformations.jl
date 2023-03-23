```@meta
CurrentModule = PBWDeformations
DocTestSetup  = quote
    using PBWDeformations
    using Oscar
end
```

# Smash products deformations

## Constructors
```@docs
DeformationMap
deform
symmetric_deformation
```

## SmashProductDeformLie struct
```@docs
SmashProductDeformLie
```

## Functions
The [`SmashProductDeformLie`](@ref) struct can be used as an argument for the following functions:
- `gens`
- `ngens`
- `change_base_ring`
