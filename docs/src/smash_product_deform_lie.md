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
- `gen`
- `gens`
- `ngens`
- `change_base_ring`

For `gen`, `gens`, and `ngens`, on can supply a symbol to choose the part of the smash product to use: `:L` for the Lie algebra, and `:V` for the module.
