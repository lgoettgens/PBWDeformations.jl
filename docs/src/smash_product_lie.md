```@meta
CurrentModule = PBWDeformations
DocTestSetup  = quote
    using PBWDeformations
    using Oscar
end
```

# Smash products

## Constructors

### General case
```@docs
smash_product
```

### Highest weight / GAP case
```@docs
```

## SmashProductLie struct
```@docs
SmashProductLie
```

## Functions
The [`SmashProductLie`](@ref) struct can be used as an argument for the following functions:
- `gen`
- `gens`
- `ngens`

For `gen`, `gens`, and `ngens`, on can supply a symbol to choose the part of the smash product to use: `:L` for the Lie algebra, and `:V` for the module.
