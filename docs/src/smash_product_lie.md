```@meta
CurrentModule = PBWDeformations
DocTestSetup  = quote
    using PBWDeformations
    using Oscar
end
```

# Smash products

## Constructors

### General / GAP case

```@docs
smash_product_lie
```

<!-- ### SO case

For the orthogonal Lie algebras `so_n` there is a different constructor, that results in the well-known basis of `so_n` given by `x_i_j` ``= E_{i,j} - E_{j,i}`` for ``1 \leq i < j \leq n``.

```@docs
``` -->

## SmashProductLie struct
```@docs
SmashProductLie
```

## Functions
The [`SmashProductLie`](@ref) struct can be used as an argument for the following functions:
- `gens`
- `ngens`
- `change_base_ring`
