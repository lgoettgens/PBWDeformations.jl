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
smash_product_deform_lie
smash_product_symmdeform_lie
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
